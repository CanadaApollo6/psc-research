import base64
import importlib.util
from pathlib import Path
import struct
import sys
import unittest
import xml.etree.ElementTree as ET
import zlib

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import qualify_ubash3a_spectra as spectra


class SpectrumQualificationTests(unittest.TestCase):
    def array(self, values, compressed=True):
        data = struct.pack("<" + "d" * len(values), *values)
        if compressed:
            data = zlib.compress(data)
        kind = "MS:1000574" if compressed else "MS:1000576"
        return ET.fromstring('<binaryDataArray xmlns="http://psi.hupo.org/ms/mzml">'
                             '<cvParam accession="MS:1000523"/><cvParam accession="' + kind + '"/>'
                             '<binary>' + base64.b64encode(data).decode() + '</binary></binaryDataArray>')

    def test_vendor_style_zlib_double_array(self):
        values = spectra.decode_array(self.array([101.123, 202.456]), 2)
        self.assertEqual(values.tolist(), [101.123, 202.456])

    def test_uncompressed_double_array(self):
        values = spectra.decode_array(self.array([1.0, 0.0], False), 2)
        self.assertEqual(values.tolist(), [1.0, 0.0])

    def test_truncated_array_is_rejected(self):
        with self.assertRaises(ValueError):
            spectra.decode_array(self.array([1.0]), 2)

    def test_nonfinite_fragment_is_rejected(self):
        with self.assertRaises(ValueError):
            spectra.decode_array(self.array([float("nan")]), 1)


if __name__ == "__main__":
    unittest.main()
