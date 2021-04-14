# -*- coding: utf-8 -*-
#BEGIN_HEADER
from MetabolicModelGapfilling.core.gapfillingmodule import GapfillingModule
#END_HEADER


class MetabolicModelGapfilling:
    '''
    Module Name:
    MetabolicModelGapfilling

    Module Description:
    A KBase module: MetabolicModelGapfilling
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.module = GapfillingModule(config)
        #END_CONSTRUCTOR
        pass


    def gapfill_metabolic_model(self, ctx, params):
        """
        :param params: instance of unspecified object
        :returns: instance of unspecified object
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN gapfill_metabolic_model
        output = self.module.run_gapfill_metabolic_model(params,ctx)
        #END gapfill_metabolic_model

        # At some point might do deeper type checking...
        if not isinstance(output, object):
            raise ValueError('Method gapfill_metabolic_model return value ' +
                             'output is not type object as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
