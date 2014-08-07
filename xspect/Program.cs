using System;
using System.Collections.Generic;
using System.Text;

namespace xspect
{
    class Program
    {
    
    	static constant string VERSION = "0.0.1";
    	static constant double PROTON = 1.007276466;
        static void Main(string[] args)
        {
            if (args.Length < 1)
            {
                showUsage();
                System.Console.Error.WriteLine("not enough arguments for xspect!");
                return;
            }
            System.Globalization.CultureInfo ci = new System.Globalization.CultureInfo("en-us");
            string rawfile_name = "";
            string peptides_filename = "";  // peptide list in TSV format for quantification
            // set output to stdout by default (can be changed to FileStreamWriter with -o parameter)
            System.IO.StreamWriter output = new System.IO.StreamWriter(Console.OpenStandardOutput());
            int q = 10;     // pick 'q' peaks per 100Da
            double min_intensity_cutoff = 0.0;  // peak picking threshold
            for(int i=0; i<args.Length; i++)
            {
                switch (args[i])
                {
                    case "-h":
                    case "--help":
                        showUsage(); return;
                    case "-o":
                    case "--output":
                        output = new System.IO.StreamWriter(args[++i]); break;
                    case "-q":
                        q = Int16.Parse(args[++i]); break;
                    case "-f":
                    case "--format":
                        ++i;
                        break;
                    case "-r":
                    case "--ratios":
                        peptides_filename = args[++i]; break;
                    default:
                        rawfile_name = args[i]; break;
                }
            }

            if (rawfile_name.Equals(""))
            {
                showUsage(); return;
            }

            
            MSFileReaderLib.IXRawfile5 raw = (MSFileReaderLib.IXRawfile5) new MSFileReaderLib.MSFileReader_XRawfile();
            raw.Open(rawfile_name);
            raw.SetCurrentController(0, 1);
            int first_scan = new int();
            raw.GetFirstSpectrumNumber(ref first_scan);
            int last_scan = new int();
            raw.GetLastSpectrumNumber(ref last_scan);
            object all_filters = null;
            int num_filters = 0;
            raw.GetFilters(ref all_filters, ref num_filters);
            string ms1_filter = null;
            string[] ms2_filters = new string[num_filters];
            for (int f = 0; f < num_filters; f++)
            {
                string filter = ((string[])all_filters)[f];
                if (filter.Contains("Full ms"))
                {
                    ms1_filter = filter;
                }
            }

            for (int i = first_scan; i <= last_scan; i++)
            {
                string filter = null;
                raw.GetFilterForScanNum(i, ref filter);
                if (filter == null || filter.Equals("") || filter.Equals(ms1_filter))
                {
                    
                }
                else
                {
                    double pre_mass = new double();
                    raw.GetPrecursorMassForScanNum(i, 2, ref pre_mass);
                    object labels = null;
                    object values = null;
                    int num_values = new int();
                    raw.GetTrailerExtraForScanNum(i, ref labels, ref values, ref num_values);
                    int charge = 0;
                    string[] label_strings = (string[]) labels;
                    string[] value_strings = (string[]) values;
                    for (int v = 0; v < num_values; v++)
                    {
                        //output.WriteLine(label_strings[v] + " " + value_strings[v]);
                        if (label_strings[v].Equals("Charge State:"))
                        {
                            charge = Int16.Parse(value_strings[v]);
                        }
                    }
                    //pre_mass = pre_mass * charge - (charge * PROTON)
                    output.WriteLine("BEGIN IONS");
                    output.WriteLine("TITLE={0} {1:D}", rawfile_name, i);
                    output.WriteLine("CHARGE={0:D}+", charge);
                    output.WriteLine("PEPMASS={0}", pre_mass.ToString(ci));
                    //output.WriteLine("ISO={0:D}", 0);	// Isotopic precursor
                    object data = null;
                    object flags = null;
                    int num_peaks = new int();
                    double centroid = new double();
                    raw.GetMassListFromScanNum(ref i, null, 0, 0, 0, 0, ref centroid, ref data, ref flags, ref num_peaks);
                    if (data == null)
                    {
                        continue;
                    }
                    //double[,] peaks = (double[,])data;
                    SortedDictionary<double, double> peaks = pickPeaks((double[,])data, q, min_intensity_cutoff);
                    foreach (KeyValuePair<double, double> p in peaks)
                    {
                        output.WriteLine("{0} {1}", p.Key.ToString(ci), p.Value.ToString("0.00",ci));
                    }
                    output.WriteLine("END IONS");
                    output.WriteLine("");
                }

            }

            raw.Close();
        }

        private static void showUsage()
        {
            System.Console.Error.Write(USAGE);
        }

        private static SortedDictionary<double, double> pickPeaks(double[,] peaks, int q, double threshold)
        {
            SortedDictionary<double, double> intens_map = new SortedDictionary<double, double>();
            SortedDictionary<double, double> peak_map = new SortedDictionary<double,double>();
            double lowerbound = 0.0;
            for (int i = 0; i < peaks.Length / 2; i++)
            {
                double mz = peaks[0, i];
                double intens = peaks[1, i];
                if (intens <= 0)
                {
                    continue;
                }
                if (mz - 100 > lowerbound)
                {
                    SortedDictionary<double, double>.Enumerator e = intens_map.GetEnumerator();
                    int count = 0;
                    while(e.MoveNext() && q>=count++) {
                        peak_map.Add(e.Current.Value, e.Current.Key);
                    }
                    intens_map = new SortedDictionary<double, double>();
                    lowerbound = mz;
                }
                while (intens_map.ContainsKey(intens))
                {
                    intens += 0.01;	// slightly favors later peaks in the spectrum
                }
                intens_map.Add(intens, mz);
            }
            SortedDictionary<double, double>.Enumerator e2 = intens_map.GetEnumerator();
            int count2 = 0;
            while (e2.MoveNext() && q * peaks[0, peaks.Length/2 - 1] - lowerbound / 100 >= count2++)
            {
                peak_map.Add(e2.Current.Value, e2.Current.Key);
            }
            return peak_map;
        }
        
    private const string USAGE = "Xspect - version "+VERSION"\n"+
    		"Extracts MS2 spectra from binary raw files of Thermo Scientific"+ 	
    		"instruments and converts them into the MGF format.\n\n"+
    	    "USAGE:\n"+
			"\txspect.exe [OPTIONS] <RAWFILE>\n"+
            "OPTIONS:\n"+
            "\t-h --help\t\tshow help\n"+
            "\t-q <NUMBER>\t\tpick 'q' most intens peaks per 100Da (DEFAULT=10)\n"+
            "\t-o --output <FILE>\tdirect output to file\n"+
			"EXAMPLE:\n"+
			"\txspect.exe -q 6 -o output.mgf rawfile.raw\n";
    }
}
