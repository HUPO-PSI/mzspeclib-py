"""Example script to write out the first n spectra from a library"""
import click

from mzspeclib import SpectrumLibrary
from mzspeclib.spectrum_library import FormatInferenceFailure
from mzspeclib.backends import TextSpectralLibraryWriter
from mzspeclib.cluster import SpectrumCluster
from mzspeclib.index import MemoryIndex, SQLIndex
from mzspeclib.spectrum import Spectrum

@click.command('first_n_entries')
@click.argument('inpath', type=click.Path(exists=True))
@click.option("-i", "--input-format", type=click.Choice(sorted(SpectrumLibrary.supported_file_extensions())),
              default=None)
@click.option("-n", '--spectra-to-read', type=int, default=20)
def main(inpath, input_format, spectra_to_read: int=20):
    """Read only the first `n` spectra from the input library and write them to STDOUT in text format"""
    if SQLIndex.exists(inpath):
        index_type = SQLIndex
    else:
        index_type = MemoryIndex
    click.echo(f"Opening {inpath}", err=True)
    try:
        library = SpectrumLibrary(filename=inpath, index_type=index_type, format=input_format, create_index=False)
    except FormatInferenceFailure as err:
        click.echo(f"{err}", err=True)
        raise click.Abort()

    stream = click.get_text_stream('stdout')
    writer = TextSpectralLibraryWriter(stream)
    writer.write_header(library)

    for i, entry in enumerate(library, 1):
        if i > spectra_to_read:
            break
        if isinstance(entry, Spectrum):
            writer.write_spectrum(entry)
        elif isinstance(entry, SpectrumCluster):
            writer.write_cluster(entry)

    writer.close()


if __name__ == "__main__":
    main.main()