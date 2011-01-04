from sage.all import db, db_save

class autodb(object):
    def __init__(self, f):
        self.f = f

    def __call__(self, *args):
        fn = ("%s" + ("_%i" * len(args))) % ((self.f.__name__, ) + args)
        try:
            return db(fn)
        except IOError:
            res = self.f(*args)
            db_save(res, fn)
            return res

