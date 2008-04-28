;;; emacs-l-helpers.el --- Set of "helping hand" functions for using emacs in
;;;                        this project

;; version=2007-10-16.19

;; Copyright (C) 2007 Stefan Bienert
;; Copyright (C) 2007 Center for Bioinformatics, University of Hamburg
;;
;; See COPYING file in the top level directory of this tree for licence.


;;; Commentary:
;;  Functions to make the work with Emacs easier...

;;; History:
;;  2007-10-10 bienert: created

;;; Code:
(defgroup elh nil "The customisation group for Emacs' lil' helper"
  :version '0.1
  :group 'Development
)

(defcustom project-tree-root nil
       "Path to the root of the projects' directory tree"
       :type 'directory
       :require 'emacs-l-helper
       :group 'elh)

(defun elh-set-project-tree-root (ptrpath)
  "Set the path to the projects' tree root.
Argument PTRPATH will be storend in project-tree-root."
  (interactive "MPath? ")
  (set-variable 'project-tree-root ptrpath)
  ;;'(project-tree-root 
  ;;  ptrpath nil (emacs-l-helper))
  )

(defun elh-check-location (pwdpath ptrpath)
  "Checks whether the path of pwdpath starts with ptrpath.
Argument pwdpath is the path to be checked.
Argument ptrpath is the path to be checked for."
  (if (numberp (string-match (expand-file-name ptrpath)
                             (expand-file-name pwdpath)))
      't
    (eval nil))
  )

(defun elh-int-put-licences_copyright_info (&optional symbol)
  "Write copyright strings and a reference to the COPYING document to buffer.
The names of copyright holders are read from file \"cpowner\", or given
interactively if not found. The years is set automatically by system time.
Optional argument SYMBOL takes a string to be writtin in front of each line."
  (interactive "MSymbol? ")
  ;; list of variables
  (defvar cpstring)
  (defvar csbuf)
  (defvar cpownerpath)
  (defvar v1)
  (when (not symbol) (setq symbol ""))
  (setq cpstring (format "%s Copyright (C) %s " symbol
                                               (format-time-string "%Y")))
  (setq csbuf (current-buffer))
  (setq cpownerpath (format "%s/cpowner" project-tree-root))
  (if (file-readable-p cpownerpath)
      (progn
        (with-temp-buffer
          (insert-file-contents cpownerpath)
          (goto-char (point-min))
          (while (not (eobp))
            (princ cpstring csbuf)
            (princ (thing-at-point 'line) csbuf)
            (forward-line))))
    (while (not (string= (
              setq v1 (read-from-minibuffer "Copyright holder (none to stop)? ")
                               ) ""))
      (princ cpstring csbuf)
      (princ v1 csbuf)
      (terpri csbuf)))
  (princ symbol csbuf)
  (terpri csbuf)
  (princ symbol csbuf)
  (princ " See COPYING file in the top level directory of this tree for licence." csbuf)
  (terpri csbuf)
  )

(defun elh-put-file-to-authors (afile)
  "Put filename into author's list in AUTHORS"
  ;; list of variables
  (defvar startp)
  (defvar ostartp)
  (defvar endp)
  (defvar fre)
  (defvar lmatch)
  (defvar currp)
  (defvar scp)
  (setq startp (point))
  (setq ostartp (point))
  (setq endp (re-search-forward "^[[:blank:]]*$" nil t))
  (when (not endp) (setq endp (point-max)))
  (goto-char endp)
  (setq fre (format " %s."(file-name-sans-extension afile)))
  (set-variable 'case-fold-search 'nil)
  (if (not (search-backward fre startp t))
      (progn
        (set-variable 'case-fold-search t)
        (setq lmatch (re-search-backward "[^[:blank:]]" startp t))
        (when (not lmatch) (setq lmatch startp))
        (when (= lmatch startp)
          ;;(terpri (current-buffer))
          (princ "       " (current-buffer))
          (setq lmatch (point))
          )
        (setq startp (line-beginning-position))
        ;;(when (>(- (+ lmatch (length afile)) startp) 78)
        ;;  (terpri (current-buffer))
        ;;  (princ "       " (current-buffer)))
        (princ " " (current-buffer))
        (princ afile (current-buffer))    
        (setq endp (re-search-forward "^[[:blank:]]*$" nil t))
        )
    (progn
      (setq fre (format "\ %s\\.\\(:?\\[[^[:blank:]]*%s[^[:blank:]]*\\]\\|%s\\)"(file-name-sans-extension afile) (file-name-extension afile) (file-name-extension afile)))
      ;;(princ (format "R: %s" fre) (current-buffer))
      ;;(goto-char ostartp)
      ;;(goto-char (- (point) 1))
      ;;(princ "I" (current-buffer))
      (when (not (re-search-forward fre endp t))
        ;; determine if already class or usual filename
        (setq fre (concat  "\ " (file-name-sans-extension afile)
                           "\\.\\[[^[:blank:]]+\\]")) ;;"\\.\\[\\w+\\]"
        (when (not (re-search-forward fre endp t))
          (setq fre (concat "\ "(file-name-sans-extension afile) "\\.\\(\\w+\\)"))
          ;; check, if search succeeded
          (if (not (re-search-forward fre endp t))
              (progn
                (error "Matching problem at AUTHORS extension.")
                (set-variable 'quit-flag t)
                )
            (progn
              ;; look out for line length
              (replace-match (concat " " (file-name-sans-extension afile)
                                     ".[\\1]") nil nil nil)
              )
            )
          )
        (goto-char (- (point) 1))
        (princ "," (current-buffer))
        (princ (file-name-extension afile) (current-buffer))
        (setq endp (re-search-forward "^[[:blank:]]*$" nil t))
        )
      ;;(princ (format "M: %s" (match-string 1)) (current-buffer))
      )
    )
  (when (not endp) (setq endp (point-max)))
  (goto-char ostartp)
  (set-variable 'fill-column 78)
  ;; sort file extensions - begin
  ;;   loop between ostaxrtp and endp
  (setq currp ostartp)
  (while (< currp endp)
    (setq scp (search-forward "[" endp 't))
    (setq currp (search-forward "]" endp 't))
   (if (or (eq scp nil) (eq currp nil))
       (setq currp endp)
     (progn
       (setq currp (- currp 1))
       (if (>= currp scp)
           (sort-regexp-fields nil "\\([^,]+\\)" "\\1" scp currp)
         (progn
           (error "Bracketstring of file extensions not in right order")
           (setq currp endp)
         )
       )
       ))
    )
  ;; sort file extensions - end
  (sort-regexp-fields nil "\ \\([^[:blank:]]+\\)" "\\1" ostartp endp)
  ;;(set-left-margin startp endp 8)
  (indent-region ostartp endp 8)
  (fill-region ostartp endp 'left nil t) ;;fill-individual-paragraphs
  (set-variable 'case-fold-search t)
  )

(defun elh-add-author-file (author lname afile)
  "Edit AUTHORS. Add new authors, add new files.
Argument symbol takes the name to be searched."
  (interactive "MAuthor? ")
  ;; list of variables
  (defvar authorpath)
  (defvar atregexp)
  (defvar atbuf)
  (defvar apoint)
  (setq authorpath (format "%s/AUTHORS" project-tree-root))
  (if (not (get-file-buffer authorpath))
      (progn
        (if (file-readable-p authorpath)
            (progn
              (setq atregexp (format "^%s[[:blank:]]+(%s):" author lname))
              (setq atbuf (find-file-noselect authorpath nil "p"))
              (save-current-buffer
                (set-buffer atbuf)
                (goto-char (point-min))
                (setq apoint (re-search-forward atregexp nil t))
                (if (not apoint)
                    (progn
                      (goto-char (point-max))
                      (terpri atbuf)
                      (princ author atbuf)(princ " (" atbuf)(princ lname atbuf)
                      (princ "):" atbuf)
                      (terpri atbuf)
                      (elh-put-file-to-authors afile)
                      )
                    (progn
                      (when (> (current-column) 0)
                        (forward-line)
                        )
                      (elh-put-file-to-authors afile)                      
                      ))
                (when (buffer-modified-p atbuf) 
                  (save-buffer 0))
                )
              (when (not (kill-buffer atbuf))
                (error "Problems at closing file AUTHORS!")))
          (error "Could not access file AUTHORS!")))
    (error "File AUTHORS is currently visited, close buffer first."))
)

(provide 'emacs-l-helpers)

;;; emacs-l-helpers.el ends here

;; Local variables:
;; eval: (add-hook 'write-file-hooks 'time-stamp)
;; time-stamp-start: ";; version="
;; time-stamp-format: "%:y-%02m-%02d.%02H"
;; time-stamp-end: "$"
;; End: