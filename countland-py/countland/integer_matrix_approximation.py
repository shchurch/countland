## INTEGER MATRIX FACTORIZATION
# 2022 SHC
# Modified from the SUSTain package: https://github.com/kperros/SUSTain
# Performs integer matrix approximation for use in reducing
#   dimensionality of single-cell seq datasets
import numpy as np


class IMA_params:
    def __init__(
        self, rank, u_bounds, l_bounds=[0, 0], maxiter=1000000.0, stop_crit=0.0001
    ):
        assert isinstance(rank, int), "rank expects an integer"

        assert isinstance(u_bounds, list), "bounds expects a list of length 2"
        assert isinstance(l_bounds, list), "bounds expects a list of length 2"
        assert len(u_bounds) == 2, "bounds expects a list of length 2"
        assert len(l_bounds) == 2, "bounds expects a list of length 2"

        ## Set parameters
        self.rank = rank

        self.u_bounds = u_bounds
        self.l_bounds = l_bounds

        self.maxiter = maxiter
        self.stop_crit = stop_crit

        self.alpha = (
            0.01  # this is only used with max inner iter, also could be removed
        )


# Function to update U, V, and Lambda each iteration
# Note: this originally used the variable V everywhere, which was confusing
#   I changed V to M
def IMA_Update_Factor(M, coeff, mkrp, mode, lambda_, params):

    r = M.shape[1]
    for k in range(r):  # for each of the ranks

        core = np.dot(M, np.multiply(lambda_, coeff[:, k]))

        ## will be subtracted to update based on new lambda[k]
        core_k = np.multiply(M[:, k], (lambda_[k] * coeff[k, k]))

        ## update scalar weight of k-th component
        delta_lambda_k = (np.dot(M[:, k].T, (mkrp[:, k] - core))) / (
            coeff[k, k] * np.power(np.linalg.norm(M[:, k]), 2)
        )
        lambda_[k] = np.max([1, np.round(lambda_[k] + delta_lambda_k)])

        ## adjust core with new lambda[k]
        core = core - core_k + (np.multiply(M[:, k], (lambda_[k] * coeff[k, k])))

        ## get k-th column of factor matrix
        Mk = M[:, k] + ((mkrp[:, k] - core) / (lambda_[k] * coeff[k, k]))

        ## project to integer constrained set
        Mk = np.round(Mk)
        if ~(np.isinf(params.l_bounds[mode])):  # unless lower bound infinite
            Mk[Mk < params.l_bounds[mode]] = params.l_bounds[
                mode
            ]  # bound values on low end
        if ~(np.isinf(params.u_bounds[mode])):  # unless upper bound infinite
            Mk[Mk > params.u_bounds[mode]] = params.u_bounds[
                mode
            ]  # bound values on high end

        ## update k-th column of factor matrix
        M[:, k] = Mk

        ## avoid zero lock
        if np.all(M[:, k] == 0):
            assert params.u_bounds[mode] >= 1 and params.l_bounds[mode] >= 0
            M[
                np.random.randint(0, M.shape[0], 1), k
            ] = 1  # add 1 to random value in k-th column

    return M, lambda_


# function to rescale if max val is above upper bound
def IMA_Compute_Init_Scaled(h, l_bound, u_bound):

    for j in range(0, h.shape[1]):
        maxval = np.max(h[:, j])

        if maxval == 0:  # if 0, stop now
            raise ValueError("Will divide by zero!")
        if maxval <= u_bound:  # if below upper bound, great
            continue

        mult_j = u_bound / maxval  # otherwise rescale by this
        h[:, j] = np.round(
            np.multiply(mult_j, h[:, j])
        )  # note subtle difference from how MATLAB rounds 0.5

    h[h < l_bound] = l_bound  # make sure low bound is low bound
    return h


# Function to initialize U, V, and Lambda
def IMA_init(X, params):

    initU = np.zeros((X.shape[0], params.rank))

    for r in range(params.rank):
        p = np.random.permutation(X.shape[0])  # randomize row numbers
        nbig = np.round((1 / params.rank) * X.shape[0]).astype(
            int
        )  # get 1/R * num cells, rounded
        initU[p[0:nbig], r] = np.random.randint(
            params.l_bounds[0], params.u_bounds[0] + 1, nbig
        )  # grab 1/R*ncells random integers between bounds
    assert (
        np.linalg.matrix_rank(initU) == params.rank
    ), "initializing U failed to achieve target rank, try again"

    rank_flag = False  # keep trying until rank == R
    max_rank_flag = 0
    while rank_flag is False:
        pat_picked = np.random.randint(
            0, X.shape[0], params.rank
        )  # random integers 1 to num cells, choose R
        initV = X[pat_picked, :].T  # pull random cells, transpose
        initV = IMA_Compute_Init_Scaled(
            initV, params.l_bounds[1], params.u_bounds[1]
        )  # rescale to within bounds

        rank_flag = np.linalg.matrix_rank(initV) == params.rank  # check that rank == R
        max_rank_flag = max_rank_flag + 1
        assert max_rank_flag < 1000, "cannot initialize V because rank is not achieved"

    initLambda = np.ones(params.rank)
    return initU, initV, initLambda


### Function to perform integer matrix factorization
def IMA(X, params):

    initU, initV, initLambda = IMA_init(X, params)

    U = initU  # initialize U
    V = initV  # initialize V
    lambda_ = initLambda  # initialize lambda

    assert sum(lambda_ >= 1) == len(
        lambda_
    )  # check that the sum of lambda >= 1 is the length of lambda

    nX = np.power(np.linalg.norm(X, "fro"), 2)  # frobenius norm

    iter_ = 0

    # convergence arrays
    e = []
    fit = []

    ## Main loop
    while iter_ <= params.maxiter:

        ## Update U and lambda
        A = np.dot(X, V)  # mkrp
        B = np.dot(V.T, V)  # coeff
        U, lambda_ = IMA_Update_Factor(U, B, A, 0, lambda_, params)

        ## Update V and lambda
        A = np.dot(X.T, U)  # mkrp
        B = np.dot(U.T, U)  # coeff
        V, lambda_ = IMA_Update_Factor(V, B, A, 1, lambda_, params)
        V_ = np.dot(V, np.diag(lambda_))

        # Check convergence
        e = np.append(
            e,
            np.sqrt(
                np.max(
                    nX
                    - (2 * np.sum(np.sum(np.multiply(V_, A))))
                    + np.sum(np.sum(np.multiply(B, (np.dot(V_.T, V_))))),
                    0,
                )
            ),
        )
        fit = np.append(fit, 1 - (e[-1] / np.sqrt(nX)))

        if iter_ > 0:
            if (
                np.absolute(fit[-2] - fit[-1]) < params.stop_crit
            ):  # check if last step showed little change, below threshold
                # print('SUSTain_M fit for each iteration:')
                # print(fit)
                break

        # iterate
        iter_ = iter_ + 1

    return U, V, lambda_
