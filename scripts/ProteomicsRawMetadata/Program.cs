using System.Text.Json;
using ThermoFisher.CommonCore.Data.Business;
using ThermoFisher.CommonCore.Data.FilterEnums;
using ThermoFisher.CommonCore.RawFileReader;

// Read vendor scan metadata once per scan. No peptide, target or spectrum score
// enters this program. The vendor assemblies remain in ignored local storage.
if (args.Length != 2)
{
    Console.Error.WriteLine("Usage: ProteomicsRawMetadata input.raw output.json");
    return 2;
}
using var raw = RawFileReaderAdapter.FileFactory(args[0]);
if (!raw.IsOpen || raw.IsError || raw.InAcquisition)
    throw new InvalidDataException("RAW is not readable as a complete acquired file");
raw.SelectInstrument(Device.MS, 1);
var counts = new SortedDictionary<int, int>();
var ms2Scans = new List<int>();
var analyzers = new SortedDictionary<string, int>();
for (int scan = raw.RunHeaderEx.FirstSpectrum; scan <= raw.RunHeaderEx.LastSpectrum; scan++)
{
    var filter = raw.GetFilterForScanNumber(scan);
    int level = (int)filter.MSOrder;
    counts[level] = counts.GetValueOrDefault(level) + 1;
    if (filter.MSOrder == MSOrderType.Ms2)
    {
        ms2Scans.Add(scan);
        string analyzer = filter.MassAnalyzer.ToString();
        analyzers[analyzer] = analyzers.GetValueOrDefault(analyzer) + 1;
    }
}
var result = new
{
    sourceFile = Path.GetFileName(args[0]),
    instrumentModel = raw.GetInstrumentData().Model,
    firstScan = raw.RunHeaderEx.FirstSpectrum,
    lastScan = raw.RunHeaderEx.LastSpectrum,
    vendorTotalSpectra = raw.RunHeaderEx.SpectraCount,
    countsByMSOrder = counts,
    ms2ScanNumbers = ms2Scans,
    ms2MassAnalyzers = analyzers,
    runtimeMinutes = raw.RunHeaderEx.EndTime - raw.RunHeaderEx.StartTime,
    metadataMethod = "Independent single-pass loop over vendor scan filters; no spectrum scoring"
};
if (counts.Values.Sum() != raw.RunHeaderEx.SpectraCount)
    throw new InvalidDataException("Vendor scan-count sum does not equal run-header total");
File.WriteAllText(args[1], JsonSerializer.Serialize(result, new JsonSerializerOptions { WriteIndented = true }) + "\n");
Console.WriteLine($"Read {counts.Values.Sum()} scans; MS2={ms2Scans.Count}; instrument={result.instrumentModel}");
return 0;
