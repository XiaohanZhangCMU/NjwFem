$machines{'r41'} = 
{
    'nbcores'    => 32,
    'mempernode' => 0,
    'execcmd'    => '',
    'template'   => 'cines.tpl',
    'script'     => 'job.par',
    'submit'     => 'llsubmit',
    'time'       => "01:00:00",
    'args'       => "",
    'argsbe'     => "",
    'bits'       => "_64bits",
    'hostarch'   => "power_ibm_aix",
    'ccprog'     => "xlc_r",
    'f77prog'    => "xlf_r",
    'mpccprog'   => "mpcc_r",
    'mcfprog'    => "mpxlf90_r",
    'ccfopt'     => "-ma -q64 -qlanglvl=extended -qarch=pwr4 -O3 -qstrict -qtune=pwr4 -s",
    'ccfdeb'     => " -ma -q64 -qlanglvl=extended -qarch=auto -g",
    'f90flags'   => "-qsuffix=cpp=f90 ",
    'ldflags'    => "-L/applications/intel/ict/3.0/cmkl/9.0/lib/64 -lmkl -L/opt/intelruntime -lifcore -lm -lrt ",
    'ccflags'    => "",
    'lkfopt'     => "-s",
    'arprog'     => "ar -X32_64 ",
    'arflags'    => "-ruv",
    'makecmd'    => "make"
}
