# GEMMA_GWAS

This KBase module is an application wrapper for the command line application _GEMMA_.

GEMMA is open source software and supported by this [repo](https://github.com/genetics-statistics/GEMMA).

GEMMA_GWAS has only a single user interface method for accepting a KBase Trait Matrix and
Variation object for running the GWAS analysis.

Currently, GEMMA_GWAS supports only output of the normal delimited GEMMA output file,
and an interactive 

## Run association analysis to discover genome location for quantified traits

The following consists of a video tutorial for the GEMMA_GWAS analysis pipeline.

### 1. Import trait and sample metadata.

Metadata object utilize the Attribute Mapping data object and the the
_Import TSV/Excel as Attribute Mapping..._ as the the import application.

[![Import metada](http://img.youtube.com/vi/vMNTVJy3uY0/0.jpg)](http://www.youtube.com/watch?v=vMNTVJy3uY0)

### 2. Import trait matrix

Using the metadata objects created in the previous step we import and validate (against the metadata files), a trait matrix
object. You will need to make sure every column header in the trait matrix exists in the trait metadata file and every
trait in the row header column exists within the trait meta object.

[![Import metada](http://img.youtube.com/vi/7WcztKadkyI/0.jpg)](http://www.youtube.com/watch?v=7WcztKadkyI)

### 3. Import variation data from staging area.

Using the Import Variation application from the VariationUtils module we can import Variant 
Call Format (VCF) files as KBase variantion objects.

_Note: Variation files are often very large and require the assistance of uploading 
and/or endpoint services to get into the KBase staging, please refer to [globus](https://www.globus.org/)
for instructions on this._

[![Import metada](http://img.youtube.com/vi/xDY4XiG5CiU/0.jpg)](http://www.youtube.com/watch?v=xDY4XiG5CiU)

### 4. Execute the GEMMA_GWAS application.

After we created the inputs for our GEMMA_GWAS application we can execute the application. You can select either a 
multivariate model or univariate model currently. The univariate model seeks to identify genomic locations associated with
each individual trait, while the mutlivariate model seeks to identify locations within the genome associated to 
each trait collectively in the supplied trait matrix object.

[![Import metada](http://img.youtube.com/vi/jB8KNqU2sFo/0.jpg)](http://www.youtube.com/watch?v=jB8KNqU2sFo)

# Setup and test

Add your KBase developer token to `test_local/test.cfg` and run the following:

```bash
$ make
$ kb-sdk test
```

After making any additional changes to this repo, run `kb-sdk test` again to verify that everything still works.

# Installation from another module

To use this code in another SDK module, call `kb-sdk install GEMMA_GWAS` in the other module's root directory.

# Help

You may find the answers to your questions in our [FAQ](https://kbase.github.io/kb_sdk_docs/references/questions_and_answers.html) or [Troubleshooting Guide](https://kbase.github.io/kb_sdk_docs/references/troubleshooting.html).
