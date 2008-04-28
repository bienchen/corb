;;; c-l-helpers.el --- skeletons for source and header files, functions and
;;;                    structs

  
;; version=2008-02-27.10
  
;; Copyright (C) 2007 Stefan Bienert
;; Copyright (C) 2007 Center for Bioinformatics, University of Hamburg

;; See COPYING file in the top level directory of this tree for licence.
  
;;; requirements

;;; Commentary:
;; 


;;; History:
;; 

;;; ToDo
;;  - Change demands on empty files: Everything with nothing or only \s and \n
;;    in it
;;  - Manage prototypes between source & header
;;  - Offer setting of library name in source/ hedaer file from current dir
;;  - module name from source file name wo extension?
;;  - Fetch projectname from somewhere
;;  - each file created (*.c, *.h) should be opened after creation
;;  - On opening check file documentation of source & headers to be consistent
;;  - automatically create/ update library header file
;;  - on open non empty file, check module/ file name to match file-name
;;  - for including generated headers in source: not <lib/header.h> but just
;;    <header.h>
;;  - saying "none" on Module name leads to "*** Module Name ***

(require 'emacs-l-helpers)
(require 'skeleton)

;;; Code:
;; C header skeleton
(defun c-header-skeleton  (&optional clh-brief clh-module clh-lib clh-project
                                     clh-author clh-short clh-editauthor)
  "Inserts a C header skeleton into current buffer. This only makes sense for
empty buffers.
Optional argument clh-brief carries brief description.
Optional argument clh-module takes a module name.
Optional argument clh-lib takes a library name.
Optional argument clh-projects takes the project name.
Optional argument clh-author takes the authors' name.
Optional argument clh-short takes the short name of the author name.
Optional argument clh-editauthor decides whether to ask to update AUTHORS or not."
  (interactive "*P")
  ;; variables - just to avoid compiler warnings AND have a nice list
  (defvar clh-hpath)
  (when (not (elh-check-location buffer-file-name project-tree-root))
    (error (format "Not in project dir %s" buffer-file-name))
    (set-variable 'quit-flag t)
    )
  (setq clh-hpath (format "%s/src" project-tree-root))
  (skeleton-proxy-new 
   '(list
     ;; start of skeleton
     "/*" ?\n
     ;; copyright information
     '(elh-int-put-licences_copyright_info " *")
     > "*/" \n
     \n
     ;; doxygen documentation stuff
     "/*" ?\n
     " ****   Documentation header   ***" \n
     > "*" \n
     > "*  @file "(file-relative-name buffer-file-name clh-hpath) \n
     > "*" \n
     '(when (not clh-brief)
        (setq clh-brief (skeleton-read "One line description? "))
        )
     > "*  @brief "clh-brief | "  *** One line description ***" \n
     > "*" \n
     '(when (not clh-module)
     (setq v1 (skeleton-read "Module name? "))
     )
     > "*  Module: "clh-module | " *** Module name ***" \n
     > "*"\n
     '(when (not clh-lib)
        (setq clh-lib (skeleton-read "Library belongings? "))
        )
     > "*  Library: "clh-lib | "*** Library belongings ***" \n
     > "*"\n
     '(when (not clh-project)
        (setq clh-project (skeleton-read "Project name? "))
        )
     > "*  Project: "clh-project | "*** Project name ***" \n
     > "*"\n
     ;; fetch login name as author and offer to user
     '(when (not clh-author)
        (setq clh-author (user-full-name))
        (setq clh-author (skeleton-read "Author? " clh-author))
     )
     > "*  @author "clh-author \n
     > "*"\n
     > "*  @date " (format-time-string "%Y-%m-%d")\n
     > "*"\n
     > "*"\n
     > "*  Revision History:"\n
     > "*         - " (format-time-string "%Y%b%d")
     '(when (not clh-short) (setq clh-short (user-login-name)))
     " "clh-short ": created"\n
     > "*"\n
     > "*/" \n
     \n
     \n
      "#ifdef __cplusplus" ?\n
      "extern \"C\" {" ?\n
      "#endif" ?\n
     ?\n
     '(setq clh-def (upcase (file-name-sans-extension 
                             (file-name-nondirectory buffer-file-name))))
     "#ifndef "clh-def"_H" ?\n
     "#define "clh-def"_H" ?\n
     \n
     \n
     ?\n
     "#endif /* "clh-def"_H */" ?\n
     ?\n
     "#ifdef __cplusplus" ?\n
     "}" ?\n
     "#endif" ?\n
     ))
  ;; finished with skeleton, now for the remianing file handling
  (when (not clh-editauthor)
    (setq clh-editauthor (y-or-n-p "Edit AUTHORS? ")))
  (when clh-editauthor
     (elh-add-author-file (user-full-name)
                          (user-login-name)
                          (file-relative-name buffer-file-name clh-hpath))
     )
)

;; create C header file
(defun cl-create-c-header (cl-file cl-brief cl-module cl-lib cl-project
                           cl-author cl-short)
  "Create a C header file from a series of arguments.
Argument cl-file denotes the name without extension.
Argument cl-brief carries one line information for the header.
Argument cl-module takes the module name.
Argument cl-lib takes the library name.
Argument cl-project takes the project name.
Argument cl-author takes the authors' name.
Argument cl-short takes the authors' short name."
  ;; list of variables
  (defvar cl-hpath)
  (defvar cl-hpathfile)
  (defvar cl-subdir)
  (defvar cl-hdrbuf)
  (setq cl-file (format "%s.h" cl-file))
  (setq cl-hpath (format "%s/src" project-tree-root))
  (setq cl-hpathfile (format "%s/%s" cl-hpath cl-file))
  (when (not (elh-check-location cl-hpathfile project-tree-root))
    (error (format "Not in project dir: %s" cl-hpathfile))
    (set-variable 'quit-flag t)
    )
  (when (file-exists-p cl-hpathfile)
    (error "File %s already exists." cl-hpathfile)
    (set-variable 'quit-flag t)    
    )
  (setq cl-subdir (file-name-directory cl-file))
  (when (and cl-subdir (not (file-exists-p(format "%s/%s" cl-hpath cl-subdir))))
    (make-directory (format "%s/%s" cl-hpath cl-subdir))
    (message (format "Created %s/%s" cl-hpath cl-subdir)))
  (setq cl-hdrbuf (find-file-noselect cl-hpathfile))
  (save-current-buffer
    (set-buffer cl-hdrbuf)
    (lock-buffer)
    (restore-buffer-modified-p 't)
    (save-buffer 0)  
    (c-header-skeleton cl-brief cl-module cl-lib cl-project cl-author cl-short
                       t)
    (if (buffer-modified-p cl-hdrbuf) 
      (save-buffer 0)
      (progn
        (error "Problem saving header file")
        (set-variable 'quit-flag t))
      )
    (unlock-buffer)
    )
  (when (not (kill-buffer cl-hdrbuf)) 
    (error "Problems at closing header file"))
  (message (format "Created %s" cl-file))
)

;; skeleton for C source files
;; - later: add source to makefile?
(define-skeleton c-source-skeleton
  "Inserts a C source skeleton into current buffer. This only makes sense for
   empty buffers."
  nil
  (when (not (elh-check-location buffer-file-name project-tree-root))
    (message (format "Not in project dir: %s" buffer-file-name))
    (set-variable 'quit-flag t)
    )
  ;; test if file AUTHORS is open
  '(when (get-file-buffer (format "%s/AUTHORS" project-tree-root))
    (error "File AUTHORS is currently visited, close buffer first.")
    (set-variable 'quit-flag t)
    )
  ;; start of skeleton
  "/*" ?\n
  ;; copyright information
  '(elh-int-put-licences_copyright_info " *")
  > "*/" \n
  \n
  ;; doxygen documentation stuff
  "/*" ?\n
  " ****   Documentation header   ***" \n
  > "*" \n
  '(setq cl-srcdir (format "%s/src" project-tree-root))
  > "*  @file "(file-relative-name buffer-file-name cl-srcdir) \n
  > "*" \n
  '(setq cl-brief (skeleton-read "One line description? "))
  '(when (equal cl-brief "") (setq cl-brief "  *** One line description ***"))
  > "*  @brief "cl-brief \n
  > "*" \n
  '(setq cl-module (skeleton-read "Module name? "))
  '(when (equal cl-module "") (setq cl-module " *** Module name ***"))
  > "*  Module: "cl-module \n
  > "*"\n
  '(setq cl-lib (skeleton-read "Library belongings? "))
  '(when (equal cl-lib "") (setq cl-lib "*** Library belongings ***"))
  > "*  Library: "cl-lib \n
  > "*"\n
  '(setq cl-project (skeleton-read "Project name? "))
  '(when (equal cl-project "") (setq cl-project "*** Project name ***"))
  > "*  Project: "cl-project \n
  '(setq v1 "")
  > "*"\n
  ;; fetch login name as author and offer to user
  '(setq cl-author (user-full-name))
  '(setq cl-author (skeleton-read "Author? " cl-author))
  > "*  @author "cl-author \n
  > "*"\n
  > "*  @date " (format-time-string "%Y-%m-%d")\n
  > "*"\n
  > "*"\n
  > "*  Revision History:"\n
  > "*         - " (format-time-string "%Y%b%d")
  '(setq cl-short (user-login-name))
  " "cl-short ": created"\n
  > "*"\n
  > "*/" \n
  \n
  \n
  ;; finished with skeleton, now for the remianing file handling
  '(when (y-or-n-p "Edit AUTHORS? ")
     (elh-add-author-file (user-full-name)
                          (user-login-name)
                          ;;(file-name-nondirectory buffer-file-name)
                          (file-relative-name buffer-file-name cl-srcdir))
     )
  '(setq cl-header "")
  '(when (y-or-n-p "Create header? ")
     (if auto-insert
         (progn
           (auto-insert-mode)
           (setq keep-autoins 't))
       (setq keep-autoins nil))
     (cl-create-c-header (file-name-sans-extension 
                             (file-relative-name buffer-file-name cl-srcdir))
                         cl-brief cl-module cl-lib cl-project cl-author cl-short
                         )
     (if keep-autoins
         (auto-insert-mode 't))
     (setq cl-header (format "#include <%s.h>"(file-name-sans-extension
                      (file-relative-name buffer-file-name cl-srcdir)))
           )
     )
  > cl-header \n
  '(goto-char (point-max))
  '(save-buffer 0)
  )

;; skeleton for C functions
;; - function to propagate class source to header
;; - keep track of prototypes in header
;; skeleton for C structs
;; modification of C:
;; - modification of file: add in authors?
;; - add history on change?

(provide 'c-l-helpers)

;;; c-l-helpers.el ends here
  
;; Local variables:
;; eval: (add-hook 'write-file-hooks 'time-stamp)
;; time-stamp-start: ";; version="
;; time-stamp-format: "%:y-%02m-%02d.%02H"
;; time-stamp-end: "$"
;; End:

