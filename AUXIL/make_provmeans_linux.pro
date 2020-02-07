pro make_provmeans_linux

  extra_cflags = '-I"/usr/local/harris/envi55/idl87/external/include"'

  MAKE_DLL, 'provmeans', 'IDL_Load', $
    compile_directory='/home/user/Documents/Projects/falcon/IDL/Extensions/mortcanty-CRCENVI-79dee6a/AUXIL', $ ; location of "provmeans.c"
    extra_cflags = extra_cflags

end