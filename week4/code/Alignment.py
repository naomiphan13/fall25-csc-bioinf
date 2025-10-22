class Alignment:
    def __init__ (self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2
        self._v = len(seq1)
        self._w = len(seq2)

    def create_score_matrix(self, gap_penalty=-2):
        v = self._v
        w = self._w
        # Initialize the scoring matrix
        if v < 0 or w < 0:
            raise ValueError ("Length cannot be nonnegative")
        elif v == 0 or w ==0:
            raise ValueError ("Sequence cannot be empty")
        else:
            score_matrix = [[0]*(w + 1) for _ in range(v + 1)]

        # Initialize the first row and column
        for i in range(v + 1):
            score_matrix[i][0] = gap_penalty * i
        for j in range(w + 1):
            score_matrix[0][j] = gap_penalty * j   
        
        return score_matrix

    def global_alignment (self, match_score=3, mismatch_score=-3, gap_penalty=-2):
        matrix = self.create_score_matrix()
        # Fill out the scoring matrix
        for i in range(1, self._v + 1):
            for j in range(1, self._w + 1):
                if self._eq_ci(self.seq1[i - 1], self.seq2[j - 1]):
                    score = match_score
                else:
                    score = mismatch_score
                matrix[i][j] = max(
                    matrix[i - 1][j] + gap_penalty,
                    matrix[i][ j - 1] + gap_penalty,
                    matrix[i - 1][j - 1] + score
                )

        # Backtracking to get the aligned sequences
        aligned_seq1 = []
        aligned_seq2 = []
        i, j = self._v, self._w

        while i > 0 or j > 0:
            current_score = matrix[i][j]

            if self._eq_ci(self.seq1[i - 1], self.seq2[j - 1]):
                score = match_score
            else:
                score = mismatch_score

            if i > 0 and j > 0 and current_score == matrix[i - 1][j - 1] + score:
                aligned_seq1.append(self.seq1[i - 1])
                aligned_seq2.append(self.seq2[j - 1])
                i -= 1
                j -= 1
            elif i > 0 and current_score == matrix[i - 1][j] + gap_penalty:
                aligned_seq1.append(self.seq1[i - 1])
                aligned_seq2.append('-')
                i -= 1
            else:
                aligned_seq1.append('-')
                aligned_seq2.append(self.seq2[j - 1])
                j -= 1

        aligned_seq1.reverse()
        aligned_seq2.reverse()

        aligned_seq1 = ''.join(aligned_seq1)
        aligned_seq2 = ''.join(aligned_seq2)

        return aligned_seq1, aligned_seq2, matrix[self._v][self._w]
        
    def local_alignment (self, match_score=3, mismatch_score=-3, gap_penalty=-2):
        alignment_matrix = [[0] * (self._w + 1) for _ in range(self._v + 1)]

        max_score = 0
        max_pos = (0, 0)

        for i in range(1, self._v + 1):
            for j in range(1, self._w + 1):
                if self._eq_ci(self.seq1[i - 1], self.seq2[j - 1]):
                    score = match_score
                else:
                    score = mismatch_score

                alignment_matrix[i][j] = max(
                    0,
                    alignment_matrix[i - 1][ j] + gap_penalty,
                    alignment_matrix[i][ j - 1] + gap_penalty,
                    alignment_matrix[i - 1][ j - 1] + score
                )

                if alignment_matrix[i][j] > max_score:
                    max_score = alignment_matrix[i][j]
                    max_pos = (i, j)

        # Backtracking to get the aligned sequences
        aligned_seq1 = []
        aligned_seq2 = []
        i, j = max_pos

        while i > 0 and j > 0 and alignment_matrix[i][j] > 0:
            current_score = alignment_matrix[i][j]

            if self._eq_ci(self.seq1[i - 1], self.seq2[j - 1]):
                score = match_score
            else:
                score = mismatch_score

            if current_score == alignment_matrix[i - 1][j - 1] + score:
                aligned_seq1.append(self.seq1[i - 1])
                aligned_seq2.append(self.seq2[j - 1])
                i -= 1
                j -= 1
            elif current_score == alignment_matrix[i - 1][j] + gap_penalty:
                aligned_seq1.append(self.seq1[i - 1])
                aligned_seq2.append('-')
                i -= 1
            else:
                aligned_seq1.append('-')
                aligned_seq2.append(self.seq2[j - 1])
                j -= 1

        aligned_seq1.reverse()
        aligned_seq2.reverse()

        aligned_seq1 = ''.join(aligned_seq1)
        aligned_seq2 = ''.join(aligned_seq2)

        return aligned_seq1, aligned_seq2, max_score

    def fitting_alignment (self, match_score=3, mismatch_score=-3, gap_penalty=-2):
        if len(self.seq1) <= len(self.seq2):
            short_seq = self.seq1
            long_seq = self.seq2
        else:
            short_seq = self.seq2
            long_seq = self.seq1

        s = len(short_seq)
        l = len(long_seq)

        # Initial three alignment matrices - lower, middle, and upper - filled with negative infinity
        if s < 0 or l < 0:
            raise ValueError ("Length cannot be negative")
        elif s == 0 or l == 0:
            raise ValueError ("Sequence cannot be empty")
        
        alignment_matrix = [[0] * (l + 1) for _ in range(s + 1)]

        # Free gap at the beginning of the short sequence
        for i in range(1, s + 1):
            alignment_matrix[i][0] = i * gap_penalty
        for j in range(1, l + 1):
            alignment_matrix[0][j] = 0

        for i in range(1, s + 1):
            for j in range(1, l + 1):
                if self._eq_ci(short_seq[i - 1], long_seq[j - 1]):
                    score = match_score
                else:
                    score = mismatch_score

                alignment_matrix[i][j] = max(
                    0,
                    alignment_matrix[i - 1][j] + gap_penalty,
                    alignment_matrix[i][ j - 1] + gap_penalty,
                    alignment_matrix[i - 1][j - 1] + score
                )

        best_j = 0
        max_score = None

        for j in range(1, l + 1):
            if alignment_matrix[s][j] > alignment_matrix[s][j - 1]:
                max_score = alignment_matrix[s][j]
                best_j = j

        # Backtracking to get the aligned sequences
        aligned_seq1 = []
        aligned_seq2 = []
        i, j = s, best_j

        while i > 0 and j > 0 and alignment_matrix[i][j] > 0:
            current_score = alignment_matrix[i][j]

            if self._eq_ci(short_seq[i - 1], long_seq[j - 1]):
                score = match_score
            else:
                score = mismatch_score

            if current_score == alignment_matrix[i - 1][j - 1] + score:
                aligned_seq1.append(short_seq[i - 1])
                aligned_seq2.append(long_seq[j - 1])
                i -= 1
                j -= 1
            elif current_score == alignment_matrix[i - 1][j] + gap_penalty:
                aligned_seq1.append(short_seq[i - 1])
                aligned_seq2.append('-')
                i -= 1
            else:
                aligned_seq1.append('-')
                aligned_seq2.append(long_seq[j - 1])
                j -= 1

        aligned_seq1.reverse()
        aligned_seq2.reverse()

        aligned_seq1 = ''.join(aligned_seq1)
        aligned_seq2 = ''.join(aligned_seq2)

        return aligned_seq1, aligned_seq2, max_score

    def affine_alignment (self, match_score=3, mismatch_score=-3, gap_open_penalty=-5, gap_extension_penalty=-1):
        v = self._v
        w = self._w
        # Initial three alignment matrices - lower, middle, and upper - filled with negative infinity
        if v < 0 or w < 0:
            raise ValueError ("Length cannot be negative")
        elif v == 0 or w == 0:
            raise ValueError ("Sequence cannot be empty")
        
        neg_inf = float("-inf")
        lower = [[neg_inf]*(w + 1) for _ in range(v + 1)]
        middle = [[neg_inf]*(w + 1) for _ in range(v + 1)]
        upper = [[neg_inf]*(w + 1) for _ in range(v + 1)]

        middle[0][0] = 0

        # Initialize the first column
        if v >= 1:
            lower[1][0] = middle[0][0] + gap_open_penalty + gap_extension_penalty
        for i in range(2, v + 1):
            lower[i][0] = lower[i - 1][0] + gap_extension_penalty

        # Initialize the first row
        if w >=1:
            upper[0][1] = middle[0][0] + gap_open_penalty + gap_extension_penalty
        for j in range(2, w + 1):
            upper[0][j] = upper[0][j - 1] + gap_extension_penalty
        
        # Fill the alignment matrix
        for i in range(1, v + 1):
            for j in range(1, w + 1):
                if self._eq_ci(self.seq1[i - 1], self.seq2[j - 1]):
                    score = match_score
                else:
                    score = mismatch_score

                lower[i][j] = max(
                    lower[i - 1][j] + gap_extension_penalty,
                    middle[i - 1][j] + gap_open_penalty + gap_extension_penalty
                )

                upper[i][j] = max(
                    upper[i][j - 1] + gap_extension_penalty,
                    middle[i][j - 1] + gap_open_penalty + gap_extension_penalty
                )

                middle[i][j] = max(
                    lower[i][j],
                    middle[i - 1][j - 1] + score,
                    upper[i][j]
                )
        
        # Backtrack to get the aligned sequences
        aligned_seq1 = []
        aligned_seq2 = []

        i, j = self._v, self._w

        # Define the current state
        if middle[i][j] >= lower[i][j] and middle[i][j] >= upper[i][j]:
            current_state = "M"
        elif lower[i][j] >= upper[i][j]:
            current_state = "L"
        else:
            current_state = "U"

        while i > 0 or j > 0:
            if current_state == "M":
                if self._eq_ci(self.seq1[i - 1], self.seq2[j - 1]):
                    score = match_score
                else:
                    score = mismatch_score
                if i > 0 and j > 0 and middle[i][j] == middle[i - 1][j - 1] + score:
                    aligned_seq1.append(self.seq1[i - 1])
                    aligned_seq2.append(self.seq2[j - 1])
                    i -= 1
                    j -= 1
                elif middle[i][j] == lower[i][j]:
                    current_state = "L"
                elif middle[i][j] == upper[i][j]:
                    current_state = "U"
                else:
                    raise RuntimeError("Backtrack: no valid predecessor from M")
                
            elif current_state == "L":
                if i > 0 and lower[i][j] == lower[i - 1][j] + gap_extension_penalty:
                    aligned_seq1.append(self.seq1[i - 1])
                    aligned_seq2.append("-")
                    i -= 1
                elif i > 0 and lower[i][j] == (middle[i - 1][j] + gap_extension_penalty + gap_open_penalty):
                    aligned_seq1.append(self.seq1[i - 1])
                    aligned_seq2.append("-")
                    i -= 1
                    current_state = "M"
                else:
                    print(lower[i][j], lower[i - 1][j], middle[i][j], middle[i - 1][j], upper[i][j])
                    raise RuntimeError("Backtrack: no valid predecessor for L")
                
            else:
                if j > 0 and upper[i][j] == upper[i][j - 1] + gap_extension_penalty:
                    aligned_seq1.append("-")
                    aligned_seq2.append(self.seq2[j - 1])
                    j -= 1
                elif j > 0 and upper[i][j] == middle[i][j - 1] + gap_extension_penalty + gap_open_penalty:
                    aligned_seq1.append("-")
                    aligned_seq2.append(self.seq2[j - 1])
                    j -= 1
                    current_state = "M"
                else:
                    raise RuntimeError("Backtrack: no valid predecessor from U")

        aligned_seq1.reverse()
        aligned_seq2.reverse()

        aligned_seq1 = ''.join(aligned_seq1)
        aligned_seq2 = ''.join(aligned_seq2)
        max_score = max(middle[v][w], lower[v][w], upper[v][w])

        return aligned_seq1, aligned_seq2, max_score
    
    def _eq_ci(self, a: str, b: str) -> bool:
        # Case-insensitive without mutating the originals
        return a.casefold() == b.casefold()


        

        # def print_matrix(name, mat, rlabels, clabels, max_dim=20):
        #     print(f"\n{name} ({len(mat)-1} x {len(mat[0])-1})")
        #     R = min(len(mat), max_dim + 1)
        #     C = min(len(mat[0]), max_dim + 1)
        #     # header
        #     header = "     " + " ".join(f"{c:>5}" for c in clabels[:C])
        #     print(header)
        #     for i in range(R):
        #         row = [f"{rlabels[i]:>3}"] + [f"{mat[i][j]:>5.0f}" if mat[i][j] != neg_inf else " -inf" for j in range(C)]
        #         print(" ".join(row))

        # rlabels = ["-"] + list(self.seq1)
        # clabels = ["-"] + list(self.seq2)

        # print_matrix("LOWER (vertical gaps in seq2)", lower, rlabels, clabels)
        # print_matrix("MIDDLE (main scores)", middle, rlabels, clabels)
        # print_matrix("UPPER (horizontal gaps in seq1)", upper, rlabels, clabels)

        # return lower, middle, upper








        