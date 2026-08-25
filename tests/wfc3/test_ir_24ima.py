import subprocess
import pytest

from ..helpers import BaseWFC3


class TestIR24IMA(BaseWFC3):
    """Test for WFC3/IR IMA reprocessing when ASN_TAB=NONE."""
    detector = 'ir'

    def test_ir_24ima(self):
        rootname = "iddh03cvq"
        input_file = f"{rootname}_ima.fits"

        # Prepare input file.
        self.get_input_file(input_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', input_file, '-v'])

        # Compare results
        outputs = [(f"{rootname}_ima_flt.fits",
                    f"{rootname}_ima_flt_ref.fits"),
                   (f"{rootname}_ima_ima.fits",
                    f"{rootname}_ima_ima_ref.fits")]
        self.compare_outputs(outputs)
