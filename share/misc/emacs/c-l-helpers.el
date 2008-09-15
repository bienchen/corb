;;; c-l-helpers.el --- skeletons for source and header files, functions and
;;;                    structs

  
;; version=2008-09-03.11
  
;; Copyright (C) 2007 Stefan Bienert
;; Copyright (C) 2007 Center for Bioinformatics, University of Hamburg

;; See COPYING file in the top level directory of this tree for licence.
  
;;; requirements

;;; Commentary:
;; 


;;; History:
;;          - 2008-08-11: changed including of generated headers from
;;                         <foo/bar.h> to "bar.h"
;;

;;; ToDo
;;  - combine duplicated parts in c-header and c-source skeleton?
;;     - check impl. of author funcs, how to check interactive mode?
;;     - next step: build switch to ask for author editing into c-file-skeleton
;;       - two vars: One for decision of c-source,
;;         one to signal interactive mode
;;       - can we determine interactive mode automatically?
;;     - next step: substitute c-file-skeleton into c-header-skeleton
;;       -> subst. c-header into c-create-header?
;;     - used by c-source and c-header
;;     - c-header gets extension for macro/ ifdef stuff
;;     - c-source gets extension to create header
;;  - Offer setting of library name in source/ header file from current dir
;;  - module name from source file name wo extension?
;;  - Fetch projectname from somewhere
;;  - each file created (*.c, *.h) should be opened after creation
;;  - on open non empty file, check module/ file name to match file-name
;;  - later: add source to makefile?
;;    - before create header
;;    - ask
;;    - keep answer for header
;;    - new question for auto-header: create... and add to library head?
;;      - "yes, create only, no"
;;    - for headers: ask add to makefile, ask to add to library header
;;    - automatically start library header
;;  - Manage prototypes between source & header
;;  - On opening check file documentation of source & headers to be consistent
;;  - adding new classes/ structs/ ... -> add to font-lock file?
;;  - elh-add-author: Repair interactive use?
;; RoadMap
;; - copyright issue:
;;   - read from file/ interactive
;;   - store in list
;;   - write immedately to give feedback for interactive
;;   - return list
;;   - use list as input for non-interactive header creation
;; - Author-stuff
;; - use new function in c-file...
;; - adopt to header creation?
;; - adopt to c-source

(require 'emacs-l-helpers)
(require 'skeleton)

;;; Code:
;; C file skeleton
(defun c-file-skeleton  (&optional cfs-brief cfs-module cfs-lib cfs-project
                                   cfs-author cfs-short cfs-editauthor)
  "Inserts a C file skeleton into current buffer. This only makes sense for
empty buffers.
Optional argument cfs-brief carries brief description.
Optional argument cfs-module takes a module name.
Optional argument cfs-lib takes a library name.
Optional argument cfs-project takes the project name.
Optional argument cfs-author takes the authors' name.
Optional argument cfs-short takes the short name of the author name.
Optional argument cfs-editauthor decides whether to ask to update AUTHORS or not."
  ;;(interactive "*P")
  ;; variables - just to avoid compiler warnings AND have a nice list
  (defvar cfs-hpath)
  (when (not (elh-check-location buffer-file-name project-tree-root))
    (error (format "Not in project dir %s" buffer-file-name))
    (set-variable 'quit-flag t)
    )
  (setq cfs-hpath (format "%s/src" project-tree-root))
  ;; test if file AUTHORS is open
  (when (get-file-buffer (format "%s/AUTHORS" project-tree-root))
    (error "File AUTHORS is currently visited, close buffer first.")
    (set-variable 'quit-flag t)
    )
  (skeleton-proxy-new 
   '(list
  "/*" ?\n
  ;; copyright information
  '(elh-int-put-licences_copyright_info " *")
  > "*/" \n
  \n
  ;; doxygen documentation stuff
  "/*" ?\n
  " ****   Documentation header   ***" \n
  > "*" \n
  > "*  @file "(file-relative-name buffer-file-name cfs-hpath) \n
  > "*" \n
  '(when (not cfs-brief)
     (setq cfs-brief (skeleton-read "One line description? "))
     )
  '(when (equal cfs-brief "") (setq cfs-brief "  *** One line description ***"))
  > "*  @brief "cfs-brief \n
  > "*" \n
  '(when (not cfs-module)
     (setq cfs-module (skeleton-read "Module name? "))
     )
  '(when (equal cfs-module "") (setq cfs-module " *** Module name ***"))
  > "*  Module: "cfs-module \n
  > "*"\n
  '(when (not cfs-lib)
     (setq cfs-lib (skeleton-read "Library belongings? "))
     )
  '(when (equal cfs-lib "") (setq cfs-lib "*** Library belongings ***"))
  > "*  Library: "cfs-lib \n
  > "*"\n
  '(when (not cfs-project)
     (setq cfs-project (skeleton-read "Project name? "))
     )
  '(when (equal cfs-project "") (setq cfs-project "*** Project name ***"))
  > "*  Project: "cfs-project \n
  > "*"\n
  ;; fetch login name as author and offer to user
  '(when (not cfs-author)
     (setq cfs-author (user-full-name))
     (setq cfs-author (skeleton-read "Author? " cfs-author))
     )
  > "*  @author "cfs-author \n
  > "*"\n
  > "*  @date " (format-time-string "%Y-%m-%d")\n
  > "*"\n
  > "*"\n
  > "*  Revision History:"\n
  > "*         - " (format-time-string "%Y%b%d")
  '(when (not cfs-short) (setq cfs-short (user-login-name)))
  " "cfs-short ": created"\n
  > "*"\n
  > "*/" \n
  \n
  \n
     ))
  (goto-char (point-max))
  (save-buffer 0)
  (elh-add-author-confirm cfs-author 
                          (user-login-name) 
                          (file-relative-name buffer-file-name cfs-hpath)
                          t
                          cfs-editauthor)
)

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
  ;; (interactive "*P")
  ;; variables - just to avoid compiler warnings AND have a nice list
  ;; (defvar clh-hpath)
  ;; (setq clh-hpath (format "%s/src" project-tree-root))
  (c-file-skeleton clh-brief clh-module clh-lib clh-project clh-author
                   clh-short clh-editauthor)
  (goto-char (point-max))
  (skeleton-proxy-new 
   '(list
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
     )
   )
  (save-buffer 0)
)

;; create C header file
(defun cl-create-c-header (cl-file cl-brief cl-module cl-lib cl-project
                           cl-author cl-short cl-add-to-authors)
  "Create a C header file from a series of arguments.
Argument cl-file denotes the name without extension.
Argument cl-brief carries one line information for the header.
Argument cl-module takes the module name.
Argument cl-lib takes the library name.
Argument cl-project takes the project name.
Argument cl-author takes the authors' name.
Argument cl-short takes the authors' short name.
Argument cl-add-to-authors decides whether to change AUTHORS automatically"
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
                       cl-add-to-authors)
    (unlock-buffer)
    )
  (when (not (kill-buffer cl-hdrbuf)) 
    (error "Problems at closing header file"))
  (message (format "Created %s" cl-file))
)

;; skeleton for C source files
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
  ;; finished with skeleton, now for the remaining file handling
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
                         't)
     (if keep-autoins
         (auto-insert-mode 't))
     (setq cl-header (format "#include \"%s.h\""(file-name-sans-extension
                      (file-name-nondirectory buffer-file-name)))
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

