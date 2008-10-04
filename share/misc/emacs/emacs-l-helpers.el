;;; emacs-l-helpers.el --- Set of "helping hand" functions for using emacs in
;;;                        this project

;; version=2008-10-04.13

;; Copyright (C) 2007 Stefan Bienert
;; 
;; This file is part of CoRB.
;; 
;; CoRB is free software: you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or
;; (at your option) any later version.
;; 
;; CoRB is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.
;; 
;; You should have received a copy of the GNU General Public License
;; along with CoRB.  If not, see <http://www.gnu.org/licenses/>.


;;; Commentary:
;;  Functions to make the work with Emacs easier...

;;; History:
;;  2007-10-10 bienert: created
;;  2008-08-13 bienert: added auto-insert advice

(require 'autoinsert)

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
  (defvar ptrpath-tmp)
  (if (string= (substring ptrpath -1) "/")
      (setq ptrpath-tmp ptrpath)
      (setq ptrpath-tmp (format "%s/" ptrpath)))
  (if (numberp (string-match (expand-file-name ptrpath-tmp)
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
  (princ " This file is part of CoRB." csbuf)
  (terpri csbuf)
  (princ symbol csbuf)
  (terpri csbuf)
  (princ symbol csbuf)
  (princ " CoRB is free software: you can redistribute it and/or modify" csbuf)
  (terpri csbuf)
  (princ symbol csbuf)
  (princ " it under the terms of the GNU General Public License as published by"
         csbuf)
  (terpri csbuf)
  (princ symbol csbuf)
  (princ " the Free Software Foundation, either version 3 of the License, or"
         csbuf)
  (terpri csbuf)
  (princ symbol csbuf)
  (princ " (at your option) any later version." csbuf)
  (terpri csbuf)
  (princ symbol csbuf)
  (terpri csbuf)
  (princ symbol csbuf)
  (princ " CoRB is distributed in the hope that it will be useful," csbuf)
  (terpri csbuf)
  (princ symbol csbuf)
  (princ " but WITHOUT ANY WARRANTY; without even the implied warranty of"
         csbuf)
  (terpri csbuf)
  (princ symbol csbuf)
  (princ " MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" csbuf)
  (terpri csbuf)
  (princ symbol csbuf)
  (princ " GNU General Public License for more details." csbuf)
  (terpri csbuf)
  (princ symbol csbuf)
  (terpri csbuf)
  (princ symbol csbuf)
  (princ " You should have received a copy of the GNU General Public License"
         csbuf)
  (terpri csbuf)
  (princ symbol csbuf)
  (princ " along with CoRB.  If not, see <http://www.gnu.org/licenses/>." csbuf)
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
  (interactive "MAuthor?\nMShort name?\nMFile to add?")
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

;; add author name to AUTHORS with question for confirmation
(defun elh-add-author-confirm (eaac-author eaac-lname eaac-afile
                               &optional eaac-do-not-ask eaac-answer)
  "First ask and then add a new file to AUTHORS. Asking can be circumvented by
setting optional argument eaac-do-not-ask. The answer is then fetched from
eaac-answer"
  (when (not eaac-do-not-ask)
    (setq eaac-answer (y-or-n-p "Edit AUTHORS? ")))
  (when eaac-answer
    (elh-add-author-file eaac-author eaac-lname eaac-afile)
     )
  eaac-answer
  )

;; add an advice to auto-insert to ignore whitespaces & newlines
;; has to be activated via "(ad-activate 'auto-insert)" in the emacs
;; configuration file
;; we use flag 'preactivate' because info page says so for byte compilation
;; we will try it with and without this flag
(defadvice auto-insert (before elh-empty-whitespace-only-buffers preactivate)
  "If a buffer only contains whitespaces (including newlines), they will be deleted so auto-insert accepts the buffer to be empty."
  (if (not (re-search-forward "[^[:blank:]\n]" nil t))
      (delete-region (point-min) (point-max))
  )
)

(provide 'emacs-l-helpers)

;;; emacs-l-helpers.el ends here

;; Local variables:
;; eval: (add-hook 'write-file-hooks 'time-stamp)
;; time-stamp-start: ";; version="
;; time-stamp-format: "%:y-%02m-%02d.%02H"
;; time-stamp-end: "$"
;; End: