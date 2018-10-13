

class CodesRouter:

    """
    A router to control all database operations on models in the
    codes application.
    """

    def db_for_read(self, model, **hints):
        """
        Attempts to read codes models go to codes_db.
        """
        if model._meta.app_label == 'codes':
            return 'codes_db'
        return None

    def db_for_write(self, model, **hints):
        """
        Attempts to write codes models go to codes_db.
        """
        if model._meta.app_label == 'codes':
            return 'codes_db'
        return None

    def allow_relation(self, obj1, obj2, **hints):
        """
        Allow relations if a model in the codes app is involved.
        """
        if obj1._meta.app_label == 'codes' or obj2._meta.app_label == 'codes':
           return True
        return None

    def allow_migrate(self, db, app_label, model_name=None, **hints):
        """
        Make sure the codes app only appears in the 'codes_db'
        database.
        """
        if app_label == 'codes':
            return db == 'codes_db'
        return None