;+
; Open an url in the default web browser. On Windows and Mac this is easy. On
; UNIX platforms, the first time this routine is called it will ask for the
; location of your preferred web browser and save this location in APP_USER_DIR
;
; @param url {in}{required}{type=string} url to goto in the default web browser
; @requires IDL 6.1
; @author Michael Galloy, 2006
;-
pro mg_open_url, url
    compile_opt strictarr

    ; launch the default web browser with the url, unfortunately, this is
    ; platform dependent
    case !version.os_family of
        'Windows' : spawn, 'start ' + url, /hide, /nowait
        else : begin
            ; Mac OS X has a nice way of doing this...
            if (!version.os_name eq 'Mac OS X') then begin
                spawn, 'Open ' + url
                return
            endif

            ; ...but the other UNIX platforms don't
            app_readme_text = $
              ['This is the configuration directory for MG_OPEN_URL ', $
               'routine. It is used to save the location of the default ', $
               'web browser between MG_OPEN_URL invocations on UNIX ', $
               'platforms.', $
               '', $
               'It is safe to remove this directory, as it', $
               'will be recreated on demand. Note that all', $
               'settings (e.g. smoke injection depth, juicitron', $
               'angle, etc.) will revert to their default settings.']

            prefdir = app_user_dir('mg', $
                               'Michael Galloy', $
                               'default-browser', $
                               'Default browser location', $
                               app_readme_text, 1)
            preffile = filepath('default-browser', root=prefdir)

            ; if
            if (file_test(preffile)) then begin
                openr, lun, preffile, /get_lun
                browser = ''
                readf, lun, browser
                free_lun, lun
                spawn, browser + ' ' + url
            endif else begin
                browser_location = dialog_pickfile()
                openw, lun, preffile, /get_lun
                printf, lun, browser_location
                free_lun, lun
                msg = ['Your browser location has been stored in:', '', $
                       '    ' + preffile, '']
                ok = dialog_message(msg, /info)
                spawn, browser_location + ' ' + url
            endelse
        end
    endcase
end
