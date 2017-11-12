import uuid

def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]

