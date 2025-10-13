import functools

class Cli:
    def __init__(self, *args):
        self.args = args

    def __call__(self, func):
        func.__cli__ = True
        func.__remark__ = self.args[0]
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)
        return wrapper
