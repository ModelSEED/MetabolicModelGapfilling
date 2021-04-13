/*
A KBase module: MetabolicModelGapfilling
*/

module MetabolicModelGapfilling {
	typedef UnspecifiedObject ReportResults;

	funcdef gapfill_metabolic_model(UnspecifiedObject params)
    returns (UnspecifiedObject output) authentication required;    
};
