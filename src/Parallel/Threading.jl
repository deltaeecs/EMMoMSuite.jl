using Base.Threads

"""
    num_threads()

Return the number of threads available.
"""
num_threads() = nthreads()

"""
    thread_id()

Return the current thread ID.
"""
thread_id() = threadid()

export num_threads, thread_id
