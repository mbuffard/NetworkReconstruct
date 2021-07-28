from scipy.stats import hypergeom

class PValues:

    def __init__(self, left_tail, right_tail, two_tail):
        self.left_tail = left_tail
        self.right_tail = right_tail
        self.two_tail = two_tail

    def __repr__(self):
        return "(left_tail=%.4g, right_tail=%.4g, two_tail=%.4g)" % \
            (self.left_tail, self.right_tail, self.two_tail)

    def pvalue(a_true, a_false, b_true, b_false):
    # Convert the a/b groups to study vs population.
        k = a_true
        n = a_false + a_true  # total in study.
        K = a_true + b_true
        N = K + a_false + b_false

        lm = max(0, n - (N - K))
        um = min(n, K)
        if lm == um:
            return PValues(1.0, 1.0, 1.0)

        epsilon = 1e-6
        cutoff = hypergeom.pmf(k, N, K, n)
        left_tail = 0
        right_tail = 0
        two_tail = 0
        for x in range(lm, um + 1):
            p = hypergeom.pmf(x, N, K, n)
            if x <= k:
                left_tail += p
            if x >= k:
                right_tail += p
            if p <= cutoff + epsilon:
                two_tail += p

        return PValues(min(left_tail, 1.0), min(right_tail, 1.0), min(two_tail, 1.0))

def pvalue_population(k, n, K, N):
    return PValues.pvalue(k, n-k, K-k, N-K-n+k)
