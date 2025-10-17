
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




(defun bks/find-edge (edge)
  "Moves point to first row description of edge label entered."
  (interactive "MEdge label: ")
  (let ((edge-label (concat edge " [label=")))
                                        ;(beginning-of-buffer)
    (search-backward "graph G {")
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
" first-char) (progn ;(message "it's a new line")
                     (forward-line)
                     (bks/move-point-to-rownum-helper rownum)))
          ((equal ":" first-char) (progn
                                    (forward-line)
                                    (bks/move-point-to-rownum-helper rownum)))
          ((equal "(" first-char) (progn
                                    ;(message "found bracket...testing number immediately after it")
                                    (forward-char 1)
                                    (let ((this-rownum (thing-at-point 'number t)))
                                      (if (> this-rownum rownum)
                                          (progn
                                            ;(message "found the right place")
                                            (backward-char)
                                            (split-line))
                                          
                                        (progn
                                          ;(message "row too small...recursing...")
                                          (forward-line)
                                          (bks/move-point-to-rownum-helper rownum)))))))))

(defun bks/find-last-node-name ()
  "Returns the name of the final node"
  (interactive)
  ; navigate to the start of the line containing the last node's name:
  (end-of-buffer)
  (search-backward "[label=")
  (beginning-of-line)
  (let ((node-name (thing-at-point 'word t)))
    (message "%s" node-name)
    node-name))

(defun bks/concat-to-string (letter-list)
  "Takes a list of individual letters; returns them as a single word"
  (setf bks/last-concat-to-string "")
  (dolist (item letter-list)
    (setf bks/last-concat-to-string (concat bks/last-concat-to-string item)))
  bks/last-concat-to-string)

(defun bks/inc-letter (letter)
  (let* ((letters '("A" "B" "C" "D" "E" "F" "G" "H" "I" "J" "K" "L" "M" "N" "O" "P" "Q" "R" "S" "T" "U" "V" "W" "X" "Y" "Z" "A"))
         (letter-pos (cl-position letter letters :test #'equal))
         (new-letter (nth (+ 1 letter-pos) letters)))
    new-letter))

(defun bks/get-next-node-name (last-node-name)
  "Returns the next node name to insert"
  (interactive)
  (let* ( ;(last-node-name (bks/find-last-node-name))
         (last-letter (car (last (split-string last-node-name "" t))))
         (next-letter (bks/inc-letter last-letter))
         (last-node-name-root (bks/concat-to-string (reverse (cdr (reverse (split-string last-node-name "" t)))))))
    (message "Last node name: %s \nLast letter: %s \nNext letter: %s \nLast node name root: %s" last-node-name last-letter next-letter last-node-name-root)
    (if (not (equal next-letter "A"))
        (progn
          (message "Not time to add another letter")
          (concat last-node-name-root next-letter))
      (progn
        (message "Time to add another letter")
        (if (equal (length (split-string last-node-name "" t)) 1)
            (progn
              (message "Returning AA")
              "AA") ; Return AA if the node name has reached Z.
          (concat (bks/get-next-node-name last-node-name-root) "A") ; otherwise, increment the root by 1 and add A to the end.
    )))))

(defun bks/make-next-node ()
  (interactive)
  (let* ((last-node-name (bks/find-last-node-name))
         (next-node-name (bks/get-next-node-name last-node-name)))
    (end-of-buffer)
    (search-backward "shape=box]")
    (end-of-line)
    (newline)
    (newline)
    (insert next-node-name " [label=\"" next-node-name ".\n\" shape=box]")
    (beginning-of-line)))

(global-set-key (kbd "C-c d") 'bks/make-next-node)


                                        ; functions for parsing in a
                                        ; row line, to try to make
                                        ; entering row info a bit
                                        ; quicker:

(defun bks/parse-rowline (rowline)
  "parses a row line into a list"
  (let* ((first-node (cadr (split-string rowline " " t)))
         (second-node (car (split-string (cadr (split-string rowline "/ " t)) " ")))
         (first-variant (string-join (cdr (split-string (cadr (split-string (car (split-string rowline " / ")) ") ")) " " t)) " "))
         (second-variant (string-join (cdr (split-string (car (cdr (split-string rowline " / "))) " ")) " "))
         (rownum (car (split-string rowline " "))))
    (list rownum first-node first-variant second-node second-variant)))

(defun bks/test-if-rowline (rowline)
  "test if string is a rowline"
  (if (equal (car (split-string rowline "" t)) "(")
      t
    nil))

(defun bks/test-if-edge-top (rowline)
  "test if string is the top of an edge block (because contains `[label=')"
  (if (cl-search "[label=" rowline)
      t
    nil))

(defun bks/test-if-edge-bottom (rowline)
  "test if string is the bottom of an edge block (because contains `shape=box]')"
  (if (cl-search "shape=box]" rowline)
      t
    nil))

(defun bks/move-to-edge-first-rowline ()
  "when point is within an edge block, move to the edge block's first rowline"
  (interactive)
  (let ((line (thing-at-point 'line t)))
    (unless (bks/test-if-edge-top line)
      (search-backward "[label="))
    (while (not (bks/test-if-rowline (thing-at-point 'line t)))
      (forward-line))))

(defun bks/insert-edge-note-from-rowline(node1 node2)
     (interactive)
     (search-backward "[label=")
     (beginning-of-line)
     (forward-line)
     (let ((rownum (read-string "Row number: ")))
       (bks/move-point-to-rownum-helper (string-to-number rownum))
       (insert "(" rownum ") " node1 " "))
     (let* ((node1prompt (concat node1 "'s variant: "))
            (reading1 (read-string node1prompt)))
       (insert reading1 " / " node2 " "))
     (let* ((node2prompt (concat node2 "'s variant: "))
            (reading2 (read-string node2prompt)))
       (insert reading2)))

(defun bks/find-or-make-a-rowline ()
  "Find the first rowline below the edge heading, or make one"
  (interactive)
  (search-backward "[label=")
  (beginning-of-line)
  (forward-line)
  (while (or (equal "
" (substring (thing-at-point 'line t) 0 1))
             (equal ":" (substring (thing-at-point 'line t) 0 1)))
    (forward-line)) ; get past the headings
  (let ((line (thing-at-point 'line t)))
    (if (bks/test-if-edge-bottom line) ; if you have now reached the bottom:
        ;; get the two node names and bung in the node-line:
        (progn
          (let ((node1 (read-string "Node 1: "))
                (node2 (read-string "Node 2: ")))
            (bks/insert-edge-note-from-rowline node1 node2)))
      ;; you should be on a rowline, but let's just test in case you're not:
      (if (bks/test-if-rowline line)
          (bks/insert-edge-note-better)
        (message "I got confused and don't know where I am."))
              )))

(defun bks/insert-edge-note-better ()
  "Better way to insert edge notes: tries to prefill the values of the current node, and prompts if not possible"
  (interactive)
  (let* ((line (thing-at-point 'line t)))
    (message "%s" line)
    ;; if you're on a rowline:
    (if (bks/test-if-rowline line)
        (let* ((rowline (bks/parse-rowline(thing-at-point 'line t)))
               (node1 (nth 1 rowline))
               (node2 (nth 3 rowline)))
          (message "rowline: %s node1: %s node2: %s" rowline node1 node2)
          (bks/insert-edge-note-from-rowline node1 node2))
      ;; otherwise, try to find a rowline:
      (bks/find-or-make-a-rowline))))

(global-set-key (kbd "C-c e") 'bks/insert-edge-note-better)
