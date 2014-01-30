from django.template import loader, Context, Template, RequestContext

def render(name, **kwds):
    tmpl = loader.get_template(name)
    cont = Context(kwds)
    page = tmpl.render(cont)
    return page

