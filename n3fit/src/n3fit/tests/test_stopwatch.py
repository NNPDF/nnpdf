""" Tests the stopwatch does what is supposed to do """

from n3fit.stopwatch import StopWatch


def time_comparer(internal_dict, computed_dict, base_time):
    for key, timing in internal_dict.items():
        diff = timing - base_time
        computed = computed_dict[key]
        assert computed == diff


def test_register_times():
    watch = StopWatch()
    watch.register_times("now")
    time_dict = watch.stop()
    # Check that the internal times are consistent
    internal_cputimes = watch._cputimes
    internal_walltimes = watch._walltimes
    start_cpu, start_wall = watch.get_times()
    time_comparer(internal_cputimes, time_dict["cputime"], start_cpu)
    time_comparer(internal_walltimes, time_dict["walltime"], start_wall)


def test_register_ref():
    base1 = "p1"
    base2 = "p2"
    watch = StopWatch()
    watch.register_times(base1)
    watch.register_ref(base2, base1)
    time_dict = watch.stop()
    # Get the absolute times of the bases we got
    b1_cpu, b1_wall = watch.get_times(base1)
    b2_cpu, b2_wall = watch.get_times(base2)
    # Compute the differences
    cpu_diff = b2_cpu - b1_cpu
    wall_diff = b2_wall - b1_wall
    # Check that the watch got the same
    keyname = f"{base1}_to_{base2}"
    assert time_dict["cputime"][keyname] == cpu_diff
    assert time_dict["walltime"][keyname] == wall_diff
