using System;
using System.Collections.Generic;
using System.Text;

namespace Test
{
    public static class RachelTest
    {
        // want to test that all proteins whose psms with base sequence !=null are being put into database
        // want to check that all modifications are localized meaning that they come from psms with full sequences !=null
        // need to make a fake MS file that goes through search and a fake DB for the PSMs
        // want a psms whose base sequence is ambiguous (=null): make sure that this does not end up in confident psms
        // want a psm whose base sequence is not ambigous but full sequence is (ptm is not localized): make sure this does not make it in DB
        // want a psm whose base sequence is not ambigous and all psms are localized: no issues should be in pruned
        // want a psm that could come from more than one protein: ensure that all options are put into the pruned DB
    }
}
