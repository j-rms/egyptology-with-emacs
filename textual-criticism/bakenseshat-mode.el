
(defun bks/insert-edge-note(edge rownum smol-node sn-reading big-node bn-reading)
     (interactive "MEdge: \nMRow: \nMSmaller node: \nMSmaller node reading: \nMGreater node: \nMGreater node reading: ")
     (when (string-match "\\`[0-9]*[1-9][0-9]*\\'" smol-node)
         (setq smol-node (concat smol-node ".")))
     (when (string-match "\\`[0-9]*[1-9][0-9]*\\'" big-node)
       (setq big-node (concat big-node ".")))
     (bks/find-edge edge)
     (bks/move-point-to-rownum-helper (string-to-number rownum))
     (insert "(" rownum ") " smol-node " " sn-reading " / " big-node " " bn-reading)
     )

(global-set-key (kbd "C-c e") 'bks/insert-edge-note)


(defun bks/find-edge (edge)
  "Moves point to first row description of edge label entered."
  (interactive "MEdge label: ")
  (let ((edge-label (concat edge " [label=")))
    (beginning-of-buffer)
    (search-forward edge-label nil nil)
    (beginning-of-line)
    (recenter-top-bottom)
    (recenter-top-bottom)
    (forward-line)))

(global-set-key (kbd "C-c M-e") 'bks/find-edge)

(defun bks/move-point-to-rownum (rownum)
  "Having moved point to first row description of edge label using `bks/find-edge', this moves point to the correct place to insert rownum"
  (interactive "NInsert row number: ")
  (bks/move-point-to-rownum-helper rownum))

(defun bks/make-row-space ()
  "makes an empty line in preparation for the new row to be entered"
  (interactive)
  (split-line))

(defun bks/move-point-to-rownum-helper (rownum)
  "Helper function for `bks/move-point-to-rownum'"
  (let ((first-char (substring (thing-at-point 'line t) 0 1)))
    (cond ((equal (point) (point-max))
           (message "reached end of buffer")) ; quit with message
                                                   ; if you have
                                                   ; reached the
                                                   ; buffer end.
          ((equal "\"" first-char) (bks/make-row-space))
          ((equal "
" first-char) (progn (message "it's a new line")
                     (forward-line)
                     (bks/move-point-to-rownum-helper rownum)))
          ((equal "(" first-char) (progn
                                    (message "found bracket...testing number immediately after it")
                                    (forward-char 1)
                                    (let ((this-rownum (thing-at-point 'number t)))
                                      (if (> this-rownum rownum)
                                          (progn
                                            (message "found the right place")
                                            (backward-char)
                                            (split-line))
                                          
                                        (progn
                                          (message "row too small...recursing...")
                                          (forward-line)
                                          (bks/move-point-to-rownum-helper rownum)))))))))
