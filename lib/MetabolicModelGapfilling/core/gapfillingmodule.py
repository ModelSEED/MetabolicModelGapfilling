# -*- coding: utf-8 -*-

from __future__ import absolute_import

from MetabolicModelGapfilling.core.basemodule import BaseModule

import cobra
import cobrakbase
import json
import csv
import logging
import optlang
import re
from optlang.symbolics import Zero, add
import cobra.util.solver as sutil
from cobrakbase.core.converters import KBaseFBAModelToCobraBuilder
from cobrakbase.Workspace.WorkspaceClient import Workspace as WorkspaceClient
from cobrakbase.core.kbase_object_factory import KBaseObjectFactory
from cobrakbase.core.fba_utilities import KBaseFBAUtilities
from cobra.core.dictlist import DictList
from cobra.core import Gene, Metabolite, Model, Reaction
from annotation_ontology_api.annotation_ontology_apiServiceClient import annotation_ontology_api

def format_direction(input):
    if input == ">":
        return "Forward"
    elif input == "<":
        return "Reverse"
    else:
        return "Reversible"

def format_gpr(input):
    proteins = []
    if "modelReactionProteins" in input:
        for protein in input["modelReactionProteins"]:
            subunits = []
            if "modelReactionProteinSubunits" in protein:
                for subunit in protein["modelReactionProteinSubunits"]:
                    genes = []
                    if "feature_refs" in subunit:
                        for gene_ref in subunit["feature_refs"]:
                            genes.append(gene_ref.split("/")[-1])
                    if len(genes) > 0:
                        if len(genes) == 1:
                            subunits.append(genes[0])
                        else:
                            subunits.append("("+" or ".join(genes)+")")
            if len(subunits) > 0:
                if len(subunits) == 1:
                    proteins.append(subunits[0])
                else:
                    proteins.append("("+" and ".join(proteins)+")")       
    if len(proteins) == 0:
        return "None"
    elif len(proteins) == 1:
        return proteins[0]
    else:
        return "("+" or ".join(proteins)+")"
    
def format_equation(input,cpdhash):
    equation = ["",""]
    for rgt in input["modelReactionReagents"]:
        cpd_id = rgt["modelcompound_ref"].split("/")[-1]
        index = 0
        if rgt["coefficient"] > 0:
            index = 1
        if len(equation[index]) > 0:
            equation[index] += " + "
        if abs(rgt["coefficient"]) != 1:
            equation[index] += "("+str(abs(rgt["coefficient"]))+") "
        equation[index] += cpdhash[cpd_id]["name"]
    return " => ".join(equation)

class GapfillingModule(BaseModule):
    def __init__(self,config):
        BaseModule.__init__(self,config,"1.0.0","MetabolicModelGapfilling")
    
    def run_gapfill_metabolic_model(self,params,ctx):
        params = self.validate_args(params,["workspace","fbamodel_id","fbamodel_output_id"],{
            "media_ids" : ["Complete"],
            "nogrow_media_ids" : [],
            "fbamodel_workspace" : params["workspace"],
            "media_workspace" : params["workspace"],
            "target_reaction" : "bio1",
            "source_fbamodel_id" : None,
            "source_fbamodel_workspace" : params["workspace"],
            "feature_ko_list" : [],
            "reaction_ko_list" : [],
            "blacklist" : [],
            "media_supplement_list" : [],
            "minimum_target_flux" : 0.1,
            "max" : 1,
            "gapfilling_annotation_sources" : [],
            "consecutive_gapfill" : 1
        })
        print("Config:"+json.dumps(self.config))
        if len(params["media_ids"]) == 0:
            params["media_ids"] = ["Complete"]
        self.initialize_call(ctx,params["workspace"],"KBaseFBA.FBAModel",params["fbamodel_output_id"])
        
        kbase_api = cobrakbase.KBaseAPI(token=ctx["token"])
        kbase_api.ws_client = self.wsclient
        modelref = params["fbamodel_workspace"]+"/"+params["fbamodel_id"]
        modelws = params["fbamodel_workspace"]
        if re.search("/",params["fbamodel_id"]):
            modelref = params["fbamodel_id"]
            modelws = None
        kbmodel = kbase_api.get_object(params["fbamodel_id"],modelws)
        model = kbase_api.get_from_ws(params["fbamodel_id"],modelws)
        src_models = []
        if params["source_fbamodel_id"] != None:
            modelws = params["source_fbamodel_workspace"]
            if re.search("/",params["source_fbamodel_id"]):
                modelws = None
            src_models.append(kbase_api.get_from_ws(params["source_fbamodel_id"],modelws)) 
        utilities = KBaseFBAUtilities(model,model,kbase_api,None,0,100,[])
        #Setting objective function
        utilities.set_objective_from_target_reaction(params["target_reaction"],params["max"])
        #Setting minimal biomass production
        utilities.convert_objective_to_constraint(params["minimum_target_flux"],None)
        #Extending model for gapfilling
        reaction_genes = self.compute_reaction_scores(kbmodel["genome_ref"],1)
        penalties = utilities.build_model_extended_for_gapfilling(1,src_models,[kbase_api.get_from_ws(kbmodel['template_ref'])],2,reaction_genes)
        #Setting minimal reaction objective
        utilities.create_minimal_reaction_objective(penalties,default_penalty = 0)
        #Running gapfill
        gfdata = {"new":{},"reversed":{}}
        media_count = 0
        if params["consecutive_gapfill"] == 1:
            for media in params["media_ids"]:
                mediaref = params["media_workspace"]+"/"+media
                media_count += 1
                if media == "Complete":
                    mediaref = "KBaseMedia/Complete"
                    utilities.apply_media_to_model(None,100,100)
                else:
                    mediaws = params["media_workspace"]
                    if re.search("/",media):
                        mediaws = None
                        mediaref = media
                    media_obj = kbase_api.get_from_ws(media,mediaws)
                    media = media_obj.info[1]
                    utilities.apply_media_to_model(media_obj,0,100)
                #Minimizing gapfilled reactions
                gapfilling_solution = model.optimize()
                gfresults = utilities.compute_gapfilled_solution(penalties)    
                #Ensuring minimality by adding binary variables
                flux_values = utilities.binary_check_gapfilling_solution(penalties,1)
                gfresults = utilities.compute_gapfilled_solution(penalties,flux_values)
                for rxn in gfresults["new"]:
                    if rxn in gfdata["new"]:
                        if gfdata["new"][rxn] != gfresults["new"][rxn] and gfdata["new"][rxn] != "=":
                            gfdata["new"][rxn] = "="
                    else:
                        gfdata["new"][rxn] = gfresults["new"][rxn]
                for rxn in gfresults["reversed"]:
                    if rxn not in gfdata["reversed"]:
                        gfdata["reversed"][rxn] = gfresults["reversed"][rxn]
                wsname = params["workspace"]
                fba_obj = self.build_fba(params["fbamodel_output_id"]+"."+media+".gf",model,flux_values,modelref,mediaref,params["target_reaction"])
                if self.config["save_output_to_kbase"] == "1":
                    kbase_api.save_object(params["fbamodel_output_id"]+"."+media+".gf", wsname,"KBaseFBA.FBA", fba_obj)
                    self.add_created_object(wsname+"/"+params["fbamodel_output_id"]+"."+media+".gf","Gapfilled FBA in "+media)
        #Adding gapfilled reactions to KBase model
        rxn_tbl = self.add_gapfilling_solution_to_kbase_model(model,kbmodel,gfdata,mediaref,reaction_genes)
        #Saving gapfilled model
        if self.config["save_output_to_kbase"] == "1":
            kbase_api.save_object(params["fbamodel_output_id"], params["workspace"],"KBaseFBA.FBAModel", kbmodel)
            self.add_created_object(wsname+"/"+params["fbamodel_output_id"],"Gapfilled model")
        gfdata = {
            'gapfilling_summary':"Model successfully gapfilled in "+str(media_count)+" media, adding "+str(len(gfdata["new"]))+" reactions and making "+str(len(gfdata["reversed"]))+" reactions reversible.",
            'reaction_tab': {
                'is_reactions': len(gfdata["new"])+len(gfdata["reversed"]) > 0,
                'reactions': json.dumps(rxn_tbl),
                'help': 'No reactions were added by the gapfilling.'
            }
        }
        self.create_report(gfdata)
        output = {}
        return self.finalize_call(output)
    
    def compute_reaction_scores(self,genome_ref,weigh_all_events_equally = 1,weights = None):
        reaction_genes = {}
        anno_api = annotation_ontology_api()
        output = anno_api.get_annotation_ontology_events({
            "input_ref" : genome_ref,
        })
        events = output["events"]
        for event in events:
            for gene in event["ontology_terms"]:
                for term in event["ontology_terms"][gene]:
                    if "modelseed_ids" in term:
                        for rxn in term["modelseed_ids"]:
                            newrxn = re.sub("^MSRXN:","",rxn)
                            if newrxn not in reaction_genes:
                                reaction_genes[newrxn] = {}
                            if gene not in reaction_genes[newrxn]:
                                reaction_genes[newrxn][gene] = 0            
                            if weigh_all_events_equally == 1 or weights == None:
                                reaction_genes[newrxn][gene] += 1
                            elif event["description"] in weights:
                                reaction_genes[newrxn][gene] += weights[event["description"]]
                            elif event["event_id"] in weights:
                                reaction_genes[newrxn][gene] += weights[event["event_id"]]
                            elif event["id"] in weights:
                                reaction_genes[newrxn][gene] += weights[event["id"]]
        return reaction_genes
    
    def build_fba(self,fba_id,model,flux_values,model_ref,media_ref,target_rxn):
        # Saving final solution as an FBA object in KBase
        fbaobj = {
            "FBABiomassVariables": [],
            "FBACompoundBounds": [],
            "FBACompoundVariables": [],
            "FBAConstraints": [],
            "FBADeletionResults": [],
            "FBAMetaboliteProductionResults": [],
            "FBAMinimalMediaResults": [],
            "FBAMinimalReactionsResults": [],
            "FBAPromResults": [],
            "FBAReactionBounds": [],
            "FBAReactionVariables": [],
            "FBATintleResults": [],
            "MFALog": "",
            "PROMKappa": 1,
            "QuantitativeOptimizationSolutions": [],
            "__VERSION__": 1,
            "additionalCpd_refs": [],
            "allReversible": 0,
            "biomassRemovals": {},
            "biomassflux_objterms": {
                "bio1": 1
            },
            "calculateReactionKnockoutSensitivity": 0,
            "comboDeletions": 0,
            "compoundflux_objterms": {},
            "decomposeReversibleDrainFlux": 0,
            "decomposeReversibleFlux": 0,
            "defaultMaxDrainFlux": 0,
            "defaultMaxFlux": 1000,
            "defaultMinDrainFlux": -1000,
            "drainfluxUseVariables": 0,
            "fbamodel_ref": model_ref,
            "findMinimalMedia": 0,
            "fluxMinimization": 1,
            "fluxUseVariables": 0,
            "fva": 0,
            "gapfillingSolutions": [],
            "geneKO_refs": [],
            "id": fba_id,
            "inputfiles": {},
            "maximizeActiveReactions": 0,
            "maximizeObjective": 1,
            "media_list_refs": [],
            "media_ref": media_ref,
            "minimizeErrorThermodynamicConstraints": 0,
            "minimize_reaction_costs": {},
            "minimize_reactions": 0,
            "noErrorThermodynamicConstraints": 0,
            "numberOfSolutions": 1,
            "objectiveConstraintFraction": 0.1,
            "objectiveValue": flux_values[target_rxn]["forward"]-flux_values[target_rxn]["reverse"],
            "other_objectives": [],
            "outputfiles": {},
            "parameters": {
                "Auxotrophy metabolite list": "",
                "Beachhead metabolite list": "",
                "minimum_target_flux": "0.01",
                "save phenotype fluxes": "0",
                "suboptimal solutions": "1"
            },
            "quantitativeOptimization": 0,
            "reactionKO_refs": [],
            "reactionflux_objterms": {},
            "simpleThermoConstraints": 0,
            "thermodynamicConstraints": 0,
            "uptakeLimits": {}
        }
        
        for varname in flux_values:
            value = flux_values[varname]["forward"]-flux_values[varname]["reverse"]
            rxn = model.reactions.get_by_id(varname)
            variable_min = rxn.lower_bound
            variable_max = rxn.upper_bound
            variable_class = 'unknown'
            #if varname in self._fva_data:
            #    variable_min, variable_max = self._fva_data[varname]
            #    variable_class = self.get_variable_class(variable_min, variable_max)
            variable_data = {
                "class": variable_class,
                "lowerBound": rxn.lower_bound,
                "max": variable_max,
                "min": variable_min,
                "upperBound": rxn.upper_bound, 
                "other_max": [],
                "other_min": [],
                "other_values": [],
                "value": value, # in kbase we assume positive uptake negative excretion
                "variableType": "flux"
            }
            variable_key = "FBAReactionVariables"
            if varname.startswith("EX_"):
                kbvar = varname[3:]
                lower = variable_data["lowerBound"]
                variable_data["lowerBound"] = -1*variable_data["upperBound"]
                variable_data["upperBound"] = -1*lower
                lower = variable_data["min"]
                variable_data["min"] = -1*variable_data["max"]
                variable_data["max"] = -1*lower
                variable_data["value"] = -1*variable_data["value"]
                variable_data["variableType"] = "drainflux"
                variable_data["modelcompound_ref"] = "~/fbamodel/modelcompounds/id/"+kbvar
                variable_key = "FBACompoundVariables"
            elif re.search("^bio\d+$",varname):
                variable_data["variableType"] = "biomassflux"
                variable_data["biomass_ref"] = "~/fbamodel/biomasses/id/"+varname            
                variable_key = "FBABiomassVariables"
            else:
                variable_data["modelreaction_ref"] = "~/fbamodel/modelreactions/id/"+varname
                variable_data["exp_state"] = 'unknown'
                variable_data["biomass_dependencies"] = []
                variable_data["coupled_reactions"] = []
                variable_data["expression"] = 0
                variable_data["scaled_exp"] = 0
            fbaobj[variable_key].append(variable_data)
        return fbaobj
    
    #Required this function to add gapfilled compounds to a KBase model for saving gapfilled model    
    def convert_cobra_compound_to_kbcompound(self,cpd,kbmodel,add_to_model = 1):
        refid = "cpd00000"
        if re.search('cpd\d+_[a-z]+',cpd.id):
            refid = cpd.id
            refid = re.sub("_[a-z]\d+$","",refid)
        cpd_data = {
            "aliases": [],
            "charge": cpd.charge,
            "compound_ref": "~/template/compounds/id/"+refid,
            "dblinks": {},
            "formula": cpd.formula,
            "id": cpd.id,
            "inchikey": "ALYNCZNDIQEVRV-UHFFFAOYSA-M",
            "modelcompartment_ref": "~/modelcompartments/id/"+cpd.id.split("_").pop(),
            "name": cpd.name,
            "numerical_attributes": {},
            "string_attributes": {}
        }
        if add_to_model == 1:
            kbmodel["modelcompounds"].append(cpd_data)
        return cpd_data

    #Required this function to add gapfilled reactions to a KBase model for saving gapfilled model    
    def convert_cobra_reaction_to_kbreaction(self,rxn,kbmodel,cpd_hash,direction = "=",add_to_model = 1,reaction_genes = None):
        rxnref = "~/template/reactions/id/rxn00000_c"
        if re.search('rxn\d+_[a-z]+',rxn.id):
            rxnref = "~/template/reactions/id/"+rxn.id
            rxnref = re.sub("\d+$","",rxnref)
        rxn_data = {
            "id": rxn.id,
            "aliases": [],
            "dblinks": {},
            "direction": direction,
            "edits": {},
            "gapfill_data": {},
            "maxforflux": 1000000,
            "maxrevflux": 1000000,
            "modelReactionProteins": [],
            "modelReactionReagents": [],
            "modelcompartment_ref": "~/modelcompartments/id/"+rxn.id.split("_").pop(),
            "name": rxn.name,
            "numerical_attributes": {},
            "probability": 0,
            "protons": 0,
            "reaction_ref": rxnref,
            "string_attributes": {}
        }
        for cpd in rxn.metabolites:
            if cpd.id not in kbmodel["modelcompounds"]:
                cpd_hash[cpd.id] = self.convert_cobra_compound_to_kbcompound(cpd,kbmodel,1)
            rxn_data["modelReactionReagents"].append({
                "coefficient" : rxn.metabolites[cpd],
                "modelcompound_ref" : "~/modelcompounds/id/"+cpd.id
            })
        if reaction_genes != None and rxn.id in reaction_genes:
            best_gene = None
            for gene in reaction_genes[rxn.id]:
                if best_gene == None or reaction_genes[rxn.id][gene] > reaction_genes[rxn.id][best_gene]:
                    best_gene = gene
            rxn_data["modelReactionProteins"] = [{"note":"Added from gapfilling","modelReactionProteinSubunits":[],"source":"Unknown"}]
            rxn_data["modelReactionProteins"][0]["modelReactionProteinSubunits"] = [{"note":"Added from gapfilling","optionalSubunit":0,"triggering":1,"feature_refs":["~/genome/features/id/"+best_gene],"role":"Unknown"}]
        if add_to_model == 1:
            kbmodel["modelreactions"].append(rxn_data)
        return rxn_data
    
    def add_gapfilling_solution_to_kbase_model(self,cobramodel,newmodel,gapfilled_reactions,media_ref,reaction_genes = None):
        gfid = None
        rxn_table = []
        if gfid == None:
            largest_index = 0
            for gapfilling in newmodel["gapfillings"]:
                current_index = int(gapfilling["id"].split(".").pop())
                if largest_index == 0 or largest_index < current_index:
                    largest_index = current_index
            largest_index += 1
            gfid = "gf."+str(largest_index)
        newmodel["gapfillings"].append({
            "gapfill_id": newmodel["id"]+"."+gfid,
            "id": gfid,
            "integrated": 1,
            "integrated_solution": "0",
            "media_ref": media_ref
        })
        cpd_hash = {}
        for cpd in newmodel["modelcompounds"]:
            cpd_hash[cpd["id"]] = cpd
        for rxn in gapfilled_reactions["new"]:
            reaction = cobramodel.reactions.get_by_id(rxn)
            kbrxn = self.convert_cobra_reaction_to_kbreaction(reaction,newmodel,cpd_hash,gapfilled_reactions["new"][rxn],1,reaction_genes)
            kbrxn["gapfill_data"][gfid] = dict()
            kbrxn["gapfill_data"][gfid]["0"] = [gapfilled_reactions["new"][rxn],1,[]]
            rxn_table.append({
                'id':kbrxn["id"],
                'name':kbrxn["name"],
                'direction':format_direction(kbrxn["direction"]),
                'gene':format_gpr(kbrxn),
                'equation':format_equation(kbrxn,cpd_hash),
                'newrxn':1
            })
        for rxn in gapfilled_reactions["reversed"]:
            for kbrxn in newmodel["modelreactions"]:
                if kbrxn["id"] == rxn:
                    kbrxn["direction"] = "="
                    rxn_table.append({
                        'id':kbrxn["id"],
                        'name':kbrxn["name"],
                        'direction':format_direction(kbrxn["direction"]),
                        'gene':format_gpr(kbrxn),
                        'equation':format_equation(kbrxn,cpd_hash),
                        'newrxn':0
                    })
                    kbrxn["gapfill_data"][gfid] = dict()
                    kbrxn["gapfill_data"][gfid]["0"] = [gapfilled_reactions["reversed"][rxn],1,[]]
        return rxn_table
    