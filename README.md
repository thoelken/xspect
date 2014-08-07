# XSpect - eXtract mass SPECTra from Thermo raw files
This programm allows the user to convert the binary raw files from Thermo Scientific mass spectrometra to easier to interpret and often more compact clear text formats.

	USAGE:
		xspect.exe [OPTIONS] <RAWFILE>
	OPTIONS:
		-h --help			show help
		-q <NUMBER>			pick 'q' most intens peaks per 100Da (DEFAULT=10)
		-o --output <FILE>	direct output to file
	EXAMPLE:
		xspect.exe -q 6 -o output.mgf rawfile.raw
		
which will only include the 6 most intens peaks per 100Da window and save the result to the `output.mgf` file. By default, the output is printed to standard output (console) and can alternatively piped to a file (e.g. `xspect.exe rawfile.raw > out.mgf`).

For now the command line tool only exports MS^2 spectra to the Mascot Generic File (MGF) format with a few slight modifications.

Output example of one entry:

	BEGIN IONS
	TITLE=ORBI-RSLC008501.RAW 8452
	CHARGE=2+
	PEPMASS=770.395935058594
	218.195343017578 22.83
	230.176162719727 729.83
	1249.244079589844 16.16
	END IONS
	
where `TITLE` contains the file name of the raw file and the scan number
`PEPMASS` reflects the monoisotopic precursor m/z value according to the instrument. Ions are listed in lines without label containing the m/z value and intensity value separated by a single whitespace. Entries are surrounded by `BEGIN IONS` and `END IONS` and seperated by empty lines. Intensity values are limited to one significant decimal place (cutoff at two decimal places, not rounded).

By default, peak picking is performed with `q=10` in 100Da.
