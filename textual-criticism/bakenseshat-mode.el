
(defun bks/insert-edge-note(rownum smol-node sn-reading big-node bn-reading)
     (interactive "MRow: \nMSmaller node: \nMSmaller node reading: \nMGreater node: \nMGreater node reading: ")
     (when (string-match "\\`[0-9]*[1-9][0-9]*\\'" smol-node)
         (setq smol-node (concat smol-node ".")))
     (when (string-match "\\`[0-9]*[1-9][0-9]*\\'" big-node)
         (setq big-node (concat big-node ".")))
     (insert "(" rownum ") " smol-node " " sn-reading " / " big-node " " bn-reading)
     )

(global-set-key (kbd "C-c e") 'bks/insert-edge-note)


(defun bks/find-edge (edge)
  "Moves point to first row description of edge label entered."
  (interactive "MEdge label: ")
  (let ((edge-label (concat edge " [label=")))
    (beginning-of-buffer)
    (search-forward edge-label nil t)
    (beginning-of-line)
    (recenter-top-bottom)
    (recenter-top-bottom)
    (forward-line)))

(global-set-key (kbd "C-c M-e") 'bks/find-edge)
