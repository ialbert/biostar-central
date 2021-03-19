'''
Spam predictor

pip install -U scikit-learn

X are the data features

y are the labels [ 0, 0, 1 ...]

Download the enron database for testing

http://www2.aueb.gr/users/ion/data/enron-spam/

Should work on any of the datasets:

seq 1 6 | parallel -j 1 wget -q -nc http://www.aueb.gr/users/ion/data/enron-spam/preprocessed/enron{}.tar.gz

'''

from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
from sklearn.pipeline import make_pipeline

from sklearn.feature_extraction.text import CountVectorizer
from sklearn.feature_extraction.text import TfidfTransformer

from sklearn.naive_bayes import GaussianNB
from sklearn.naive_bayes import MultinomialNB
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier

from sklearn.base import TransformerMixin
from itertools import *
import sys

from joblib import dump, load

MAX_WORDS = 2000


def load_model(fname="clf.joblib"):
    nb = load(fname)
    return nb


def build_model(X, y, fname="clf.joblib"):
    nb = make_pipeline(
        CountVectorizer(max_features=MAX_WORDS),
        MultinomialNB()
    )
    nb.fit(X, y)
    if fname:
        dump(nb, fname)
    return nb


def evaluate_model(X, y):
    X_train, X_test, y_train, y_test = train_test_split(X, y)
    nb = build_model(X_train, y_train)
    y_pred = nb.predict(X_test)
    rep = classification_report(y_test, y_pred)
    print(rep)


def parse_file(stream):
    for line in stream:
        for word in line.split():
            word = word.decode("utf-8", errors="ignore").strip()
            if len(word) > 3:
                yield word


def test(fname=None):
    import tarfile

    # Take filename from command line
    fname = fname or sys.argv[1]
    tar = tarfile.open(name=fname, mode='r:gz', fileobj=None)
    elems = filter(lambda t: t.isreg(), tar)
    X, y = [], []
    for info in elems:
        # Fill in the labels
        y.append(int("spam" in info.name))
        stream = tar.extractfile(info)
        content = stream.read().decode("utf-8", errors="ignore")
        X.append(content)

    spam_count = sum(filter(None, y))
    ham_count = len(y) - spam_count
    print(f"Posts: {ham_count} ham, {spam_count} spam")

    evaluate_model(X, y)


def main():
    test()


if __name__ == '__main__':
    main()
