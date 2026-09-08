/* Publication-path compatibility shared by local processes. */

def runtimeOutdir() {
    def value = java.lang.System.getProperty('tresflow.resolvedOutdir')
    if( !value ) {
        throw new IllegalStateException('TrESFlow resolved output directory is not initialized')
    }
    return value
}
