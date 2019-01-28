/*
A KBase module: GEMMA_GWAS
*/

module GEMMA_GWAS {
    /* An X/Y/Z style reference */
    typedef string obj_ref;

	typedef structure {
	    string output_ws;
		obj_ref Variation;
		obj_ref Associations;
		/*
		    TODO:
		    GWAS parameters for end-user customization of GWAS analysis.
		*/
	} GemmaGwasInput;

    typedef structure {
        string report_name;
        string report_ref;
        /*
            TODO:
            Output viewer for GWAS results
        */
    } GwasResults;


    funcdef run_GEMMA_GWAS(GemmaGwasInput params) returns (GwasResults output) authentication required;
};
