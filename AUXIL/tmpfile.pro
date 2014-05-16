; NAME:
;        TMPFILE
;
; PURPOSE:
;        Get the filename and path of a temporary file.
;
; CATEGORY:
;        IO.
;
; CALLING SEQUENCE:
;
;        Result = TMPFILE( [Prefix, Suffix, Ndigits] )
;
; OPTIONAL INPUTS:
;
;        Prefix:   The prefix of the filename, ('tmp'=default).
;
;        Suffix:   The suffix of the filename, (''=default).
;
;        Ndigits:  The number of digits to append to the Prefix, (3=default).
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;
;        SEED:     The seed used by the random number generator to determine
;             the number to append to the Prefix.
;
; OUTPUTS:
;        The file specification returned is:
;
;             tmp_dir+Prefix+NNN+Suffix
;
;        where NNN is a random LONG integer with number of digits = Ndigits.
;        This routine also searches the tmp_dir for this filename. If a
;        duplicate filename is found, then another random number is chosen.
;
;        tmp_dir is the path to the /tmp directory (see FILEPATH).
;
;        If a Suffix is provided, then a '.' is appended.
;
;
; EXAMPLE:
;        Let's create a temp filename myNNNN.ps:
;
;             tmp_file = TMPFILE('my','ps',4)
;
; MODIFICATION HISTORY:
;        Written by:    Han Wen, November 1994.
;        08-FEB-1995    Fixed minor bug, random integer -&gt; random LONG integer
;        12-JUN-1995    Check existence of TMP directory and create it if no
;                       such directory exists.
;        01-JUL-1995    Use DIR_EXIST to check existence of TMP directory.
;        07-AUG-1996    Use TMPPATH (hacked version of FILEPATH) to determine
;                       TMP directory.
;-
function TMPPATH, FILENAME, ROOT_DIR=root_dir

         ON_ERROR,2                    ;Return to caller if an error
                                       ;occurs
         path = ''

         VERSION_OS = STRLOWCASE(STRMID(!VERSION.OS,0,3))
         CASE VERSION_OS OF
           'vms': root_dir = 'SYS$LOGIN:'
           'win': begin
                   root_dir = GETENV('TMP')
                   if (root_dir eq '') then root_dir = GETENV('TEMP')
                   if (root_dir eq '') then root_dir = '\tmp'
                  end
           'mac': begin
                   root_dir = !DIR
                   if (n_params() EQ 0) then filename = "IDL Temp File"
                  end
           ELSE: begin
                   root_dir = GETENV('TMP')
                   if (root_dir eq '') then root_dir = GETENV('TEMP')
                   if (root_dir eq '') then root_dir = '/tmp'
                  end
         ENDCASE

         CASE VERSION_OS OF
           'vms': BEGIN
               IF (NOT do_tmp) THEN BEGIN
                 IF (path EQ '') THEN path = '000000'
                 path = '[' + path + ']'

               ; check for a ".]" at the end of our root directory

                 IF(( STRMID(root_dir, STRLEN(root_dir)-2, 2) ne ".]") and    $
                    ( STRMID(root_dir, STRLEN(root_dir)-1, 1) eq "]") )then   $
                    root_dir = STRMID(root_dir,0,STRLEN(root_dir)-1) +'.]'
               ENDIF
             END
           'win': BEGIN
               path = '\' + path
               IF (path NE '\') THEN path = path + '\'
             END
           'mac': BEGIN
               ; make sure the root dir ends with a separator
               IF (STRMID(root_dir, STRLEN(root_dir) - 1, 1) NE ':') THEN $
             root_dir = root_dir + ':'
             END
           ELSE: BEGIN
               path = '/' + path
               IF (path NE '/') THEN path = path + '/'
             END
         ENDCASE
         RETURN, root_dir + path + filename

END

function TMPFILE, Prefix, Suffix, Ndigits, SEED=Seed

         common TMPCOM, tmpseed

         NP = N_PARAMS()
         case NP of
              0    : begin
                        name      = 'tmp'
                        ext       = ''
                        Ndigits   = 3
                     end
              1    : begin
                        name      = Prefix
                        ext       = ''
                        Ndigits   = 3
                     end
              2    : begin
                        name      = Prefix
                        ext       = '.'+Suffix
                        Ndigits   = 3
                     end
              3    : begin
                        name      = Prefix
                        ext       = '.'+Suffix
                     end
              else : message,'Invalid parameter, Ndigits='+$
                             string(Ndigits)
         endcase

         if KEYWORD_SET( SEED ) then tmpseed = SEED
         VERSION_OS = STRLOWCASE(STRMID(!VERSION.OS,0,3))

NEWRND:  numstr    = long( 10.^Ndigits * randomu(tmpseed) )
         numstr    = STRTRIM( numstr,2 )
         while (STRLEN(numstr) lt Ndigits) do $
            numstr = '0'+numstr

         filename  = name+numstr+ext
         fpath     = TMPPATH( filename )

         existfile = FINDFILE( fpath, COUNT=nexist )
         if nexist gt 0 then goto, NEWRND

         return, fpath
end