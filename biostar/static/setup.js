/*
This script should be called before other javascript.
Contains setup and utility functions.
 */

/*
Makes string substitutions simpler:

var text = "{0}, {1}! {0}!".format("Hello", "World")

Prints:

    Hello World! Hello!

 */
String.prototype.format = String.prototype.f = function () {
    var s = this,
        i = arguments.length;

    while (i--) {
        s = s.replace(new RegExp('\\{' + i + '\\}', 'gm'), arguments[i]);
    }
    return s;
};


// Get cookie from headers.
function getCookie(name) {
    var cookieValue = null;
    if (document.cookie && document.cookie != '') {
        var cookies = document.cookie.split(';');
        for (var i = 0; i < cookies.length; i++) {
            var cookie = cookies[i].trim();
            // Does this cookie string begin with the name we want?
            if (cookie.substring(0, name.length + 1) == (name + '=')) {
                cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
                break;
            }
        }
    }
    return cookieValue;
}

// Token set by Django.
csrftoken = getCookie('csrftoken');

function csrfSafeMethod(method) {
    // These HTTP methods do not require CSRF protection.
    return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
}

// PDxmRmY7FKgIzBoVGZCXG3QDOwIMrd0J1G3

// Set up Ajax functions.
$.ajaxSetup({
    crossDomain: false, // Obviates need for sameOrigin test.
    beforeSend: function (xhr, settings) {
        if (!csrfSafeMethod(settings.type)) {
            xhr.setRequestHeader("X-CSRFToken", csrftoken);
        }
    }
});


function show_message(elem, text, status) {
    var node = $("<div class='ui basic segment notify'>{0}</div>".format(text));
    node.addClass("ui {0} message".format(status));
    elem.closest("form").addClass(status).prepend(node)
}

function popover_message(elem, message, cls, timeout) {


    timeout = typeof timeout !== 'undefined' ? timeout : 1000;
    // Only works over text area
    elem = elem.find("textarea")
    var text = $('<div class="popover"></div>');
    var tag = $(text).insertBefore(elem)
    tag.addClass(cls)
    tag.text(message)
    tag.delay(timeout).fadeOut(500, function () {
        $(this).remove()
    });
}


function popup_message(elem, msg, cls, timeout) {
    timeout = typeof timeout !== 'undefined' ? timeout : 1000;
    var text = '<div></div>'
    var tag = $(text).css('z-index', '1000').insertBefore(elem);
    tag.addClass('popover ' + cls);
    tag.text(msg);
    tag.delay(timeout).fadeOut(500, function () {
        $(this).remove()
    });
}

// Triggered on network errors.
function error_message(elem, xhr, status, text) {
    popup_message(elem, "Error! readyState=" + xhr.readyState + " status=" + status + " text=" + text, "error", timeout = 5000)
}
