from django.contrib import messages
from django.shortcuts import render, redirect
import mistune

import os
import glob


def join(*args):
    return os.path.abspath(os.path.join(*args))




def docs(request, name, docs_dir):
    patt = join(docs_dir, name) + ".*"
    files = glob.glob(patt)
    if not files:
        msg = f"Cannot be find the requested page: {name} "
        messages.error(request, msg)
        return redirect("index")
    if len(files) > 1:
        msg = f"Multiple files match: {{name}}"
        messages.warning(request, msg)
    target = files[0]
    content = open(target).read()

    # Render markdown into HTML.
    if target.endswith(".md"):
        content = mistune.markdown(content)

    title = name.replace("-", " ").replace("_", " ").title()
    context = dict(content=content, title=title)
    return render(request, 'info/doc_base.html', context=context)
