include { runtimeShellExports; shellQuoteForFastqCleanup } from '../runtime_support/main'

process VALIDATE_PREPROCESSED_SPLIT {
    tag "${meta.modality}.${splitId}"
    label 'process_single'

    input:
    tuple val(splitId), val(meta), val(spec), path(r1, stageAs: 'r1/*'), path(r2, stageAs: 'r2/*'), path(rg, stageAs: 'rg/*')

    output:
    tuple val(splitId), val(meta), path('validated'), emit: validated

    script:
    def options = groovy.json.JsonOutput.toJson(spec).bytes.encodeBase64().toString()
    def rgArgs = meta.modality == 'dna' ? "--output-rg validated/SAM_RG_Header_${splitId}.tsv" : ''
    if( rg ) rgArgs += " --input-rg ${shellQuoteForFastqCleanup(rg)}"
    """
    ${runtimeShellExports(meta)}
    mkdir validated
    "\$PYTHON3_BIN" "${projectDir}/bin/preprocessed_fastqs.py" validate \
        "${options}" ${shellQuoteForFastqCleanup(r1)} ${shellQuoteForFastqCleanup(r2)} ${rgArgs}
    mv validated.txt validated/
    """
}
