def run_mh(func, stepper, rng, start, nstep):
    """
    Run a Metropolis-Hastings sampler

    Parameters
    ----------
    func: callable
        A function to calculate the log(probability) for a set of parameters
    stepper: callable
        A function to create a new set of parameters from an input set
    rng: randome number generator
        E.g. np.random.default_rng or np.random.RandomState
    start: array or sequence
        Inital set of parameters representing a state of the system
    nstep: int
        Number of steps to run in the chain

    Returns
    -------
    trials, loglike, accepted

    trials: array
        Array of shape (ntrials, npars)
    loglike: array
        Array of shape (ntrials, )
    accepted: array
        Array of shape (ntrials, ) and type bool, set to True for
        accepted points in the chain
    """
    import numpy as np

    trials, loglike, accepted = _init_data(func=func, start=start, nstep=nstep)

    for i in range(1, nstep):
        old_state = trials[i - 1]
        old_like = loglike[i - 1]

        new_state = stepper(old_state)
        new_like = func(new_state)

        loglike_ratio = new_like - old_like

        logr = np.log(rng.uniform())

        # we take np.inf or np.nan as a sign we are out of bounds
        # and thus do not accept the point in the chain
        if np.isfinite(new_like) and (
            (new_like > old_like) | (logr < loglike_ratio)
        ):
            accepted[i] = True
            loglike[i] = new_like
            trials[i, :] = new_state

        else:
            accepted[i] = False
            loglike[i] = old_like
            trials[i, :] = old_state

    return {
        'trials': trials,
        'loglike': loglike,
        'accepted': accepted,
    }


def _init_data(func, start, nstep):
    """
    Set the trials and accept array.
    """
    import numpy as np

    start = np.array(start, dtype="f8", copy=False)
    start_loglike = func(start)
    if not np.isfinite(start_loglike):
        raise RuntimeError(
            f'starting point {start} gives log prob {start_loglike}'
        )

    npars = start.size

    trials = np.zeros((nstep, npars))
    loglike = np.zeros(nstep)
    accepted = np.zeros(nstep, dtype=bool)

    trials[0, :] = start
    loglike[0] = start_loglike
    accepted[0] = True

    return trials, loglike, accepted
