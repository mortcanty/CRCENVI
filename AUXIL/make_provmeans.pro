pro make_provmeans

 extra_cflags = '/I"c:\programme\Microsoft Visual Studio\Vc98\include"'
 extra_lflags = '/Libpath:"c:\programme\Microsoft Visual Studio\Vc98\lib"'

MAKE_DLL, 'provmeans', 'IDL_Load', compile_directory='d:\idl\projects\auxil', $
 extra_cflags = extra_cflags, extra_lflags = extra_lflags

end