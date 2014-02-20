import random, hashlib, json, base64, hmac

def make_uuid():
    "Returns a unique id"
    x = random.getrandbits(256)
    u = hashlib.md5(str(x)).hexdigest()
    return u

def encode(data, key):
    text = json.dumps(data)
    text = base64.urlsafe_b64encode(text)
    digest = hmac.new(key, text).hexdigest()
    return text, digest

def decode(text, digest, key):
    if digest != hmac.new(key, text).hexdigest():
        raise Exception("message does not match the digest")
    text = base64.urlsafe_b64decode(text)
    data = json.loads(text)
    return data