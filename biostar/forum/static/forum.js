// As seen in: https://stackoverflow.com/questions/1038746/equivalent-of-string-format-in-jquery
// Usage: 'Hello {0} and {1}, and {0}'.f("World", "John");

String.prototype.format = String.prototype.f = function () {
    var s = this,
        i = arguments.length;

    while (i--) {
        s = s.replace(new RegExp('\\{' + i + '\\}', 'gm'), arguments[i]);
    }
    return s;
};

// using jQuery
function getCookie(name) {
    var cookieValue = null;
    if (document.cookie && document.cookie != '') {
        var cookies = document.cookie.split(';');
        for (var i = 0; i < cookies.length; i++) {
            var cookie = jQuery.trim(cookies[i]);
            // Does this cookie string begin with the name we want?
            if (cookie.substring(0, name.length + 1) == (name + '=')) {
                cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
                break;
            }
        }
    }
    return cookieValue;
}

csrftoken = getCookie('csrftoken');

function csrfSafeMethod(method) {
    // these HTTP methods do not require CSRF protection
    return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
}

$.ajaxSetup({
    crossDomain: false, // obviates need for sameOrigin test
    beforeSend: function (xhr, settings) {
        if (!csrfSafeMethod(settings.type)) {
            xhr.setRequestHeader("X-CSRFToken", csrftoken);
        }
    }
});


function popup_message(elem, msg, cls, timeout) {
    timeout = typeof timeout !== 'undefined' ? timeout : 1000;
    var text = '<div></div>'
    var tag = $(text).insertBefore(elem)
    tag.addClass('popover ' + cls)
    tag.text(msg)
    tag.delay(timeout).fadeOut(500, function () {
        $(this).remove()
    });
}

// Triggered on network errors.
function error_message(elem, xhr, status, text) {
    popup_message(elem, "Error! readyState=" + xhr.readyState + " status=" + status + " text=" + text, "error", timeout = 5000)
}

function apply_vote(elem, post_uid, vote_type) {

    // Toggle the button to provide feedback
    elem.toggleClass("on");

    // The vote handler.
    vote_url = "/ajax/vote/";

    $.ajax(vote_url, {
        type: 'POST',
        dataType: 'json',
        ContentType: 'application/json',
        data: {
            'post_uid': post_uid,
            'vote_type': vote_type,
        },

        success: function (data) {
            if (data.status === 'error') {
                // Untoggle the button if there was an error
                elem.toggleClass("on")
                popup_message(elem, data.msg, data.status);
            } else {
                // Success
                //popup_message(elem, data.msg, data.status);
                // Increment the post score counter
                var score = $("#score-" + post_uid)
                var value = (parseInt(score.text()) || 0) + parseInt(data.change) || 0;
                score.text(value)
            }

        },
        error: function (xhr, status, text) {
            elem.toggleClass("on")
            error_message(elem, xhr, status, text)
        }
    });
}

function remove_trigger() {
    // Makes site messages dissapear.
    $('.remove').delay(2000).slideUp(800, function () {
        $(this).remove();
    });
}



function cancel_inplace() {

    var inplace_content = $('#new-edit');
    //var inplace_title = $('inplace-title[data-value="'+ uid +'"]');
    $('.editing-drag-off').attr('draggable', true);
    //var title = $('.editable-title[data-value="'+ uid +'"]');
    var content = $('.editable');
    //var content = $('#wmd-input-' + uid);
    var hidden = $('.hide-on-edit');
    $('.hide-on-comment').show();

    //Delete the form
    inplace_content.remove();
    //inplace_title.html("");
    // Hide the container
    // Show original content
    content.show();
    //Show any blocked element
    hidden.show();

}


function highlight_search(target, content_elem, stop_list) {

    // Find the target in the content.

    var target_list = target.replace(/\s{2,}/g, ' ').split(" ");

    $.each(target_list, function (index, value) {

        var html = content_elem.html();
        var insert = "<span class='search-highlight'>" + value + "</span>";

        var new_html = html.replace(new RegExp(value, "ig"), insert);

        // Don't highlight any stop words.
        if ($.inArray(value, stop_list) === -1 && value.length >= 3) {
            console.log(value, insert, html);
            content_elem.html(new_html);
        }

    });
}

function inplace_post_edit(elem) {

    var uid = elem.data("value");
    var hidden = $('.hide-on-edit[data-value="' + uid + '"]');
    var form_container = $('inplace[data-value="' + uid + '"]');
    var dim_elem = $('.dimm-on-edit[data-value="' + uid + '"]');
    var url = '/inplace/form/';

    // Check if other posts are being edited.
    var editing = $("#new-edit");
    $('#new-comment').remove();
    $('#add-answer').html('');
    dim_elem.dimmer('show');

    if (editing.length) {
        // Remove exiting edits
        $('.hide-on-edit').show();
        $('.editable').show();
        editing.remove();
    } else {
        // Create a new edit
        editing = $('<div id="new-edit"></div>')
    }
    form_container.after(editing);
    editing.html('');

    $.ajax(url,
        {
            type: 'GET',
            dataType: 'json',
            ContentType: 'application/json',
            data: {
                'uid': uid,
            },
            success: function (data) {
                if (data.status === 'error') {
                    popup_message(elem, data.msg, data.status, 3000);
                    dim_elem.dimmer('hide');
                } else {
                    $('.editing-drag-off[id="' + uid + '"]').attr('draggable', false);
                    //alert($('.editing-drag-off[data-value="'+ uid + '"]').attr('draggable'));
                    elem.hide();
                    hidden.hide();
                    dim_elem.dimmer('hide');
                    editing.html(data.inplace_form);
                    editing.show().find('textarea').focus();
                    var preview = $('#preview');
                    preview.find('pre').addClass(' language-python');
                    preview.find('code').addClass('  language-bash language-python');

                }
            },
            error: function (xhr, status, text) {
                error_message(elem, xhr, status, text)
            }
        })
}


function highlight(text, preview) {
    var con = markdownit({
        // ESCAPES when html=true
        html: false,
        highlight: function (str, lang) {
            if (lang && hljs.getLanguage(lang)) {
                try {
                    return '<pre class="hljs"><code>' +
                        hljs.highlight(lang, str, true).value +
                        '</code></pre>';
                } catch (__) {
                }
            }
            return '<pre class="hljs"><code>' + con.utils.escapeHtml(str) + '</code></pre>';
        }
    });

    var res = con.render(text);
    preview.html(res);
    preview.find('pre').addClass('language-bash');
    preview.find('code').addClass('language-bash');

    //Prism.highlightAll()
}

function allowDrop(ev) {
    ev.preventDefault();
    //alert(ev);
}

var children = '';
var dragged_over = '';

function drag(ev, elem_id) {
    //ev.preventDefault();
    ev.dataTransfer.setData("text", elem_id);
    //let elem = $('.droptarget[data-value="' + elem_id + '"]');

    //children = '{0}'.format(elem.data('thread'));

    ev.dataTransfer.setData("text", elem_id);
    ev.target.style.opacity = "0.4";

    //alert(elem_id);

}

function drag_leave(ev, elem) {

    //let elem = $('.droptarget[data-value="' + pid + '"]');
    elem.css('border', 'none');

    //elem.css('padding', '0');
    elem.css('opacity', '1');
    //elem.css('backgroundColor', 'none');
    dragged_over = '';
    //console.log(elem.data('value'), dragged_over);

}


function reset_drag(source, target) {
    console.log(target, source);
    source.css('opacity', '1');
    source.css('border', 'none');
    target.css('opacity', '1');
    target.css('border', 'none');

}

function drag_over(ev, elem) {
    dragged_over = '';
    //console.log(elem.attr('class'), elem.data('value'));
    dragged_over = elem.data('value');
    elem.css('border', '#c2ffc2 dashed 5px');

}


function drop(ev, elem_id) {

    ev.preventDefault();

    var source = ev.dataTransfer.getData("text");
    let source_elem = $('#' + source);
    source_elem.css('opacity', '1');
    let elem = $('.droptarget[data-value="' + dragged_over + '"]');

    //alert(dragged_over);
    $.ajax('/drag/and/drop/',
        {
            type: 'POST',
            dataType: 'json',
            ContentType: 'application/json',
            data: {
                'uid': source,
                'parent': dragged_over,

            },
            success: function (data) {
                if (data.status === 'error') {
                    popup_message(elem, data.msg, data.status, 2000);
                    reset_drag(elem, source)
                } else {
                    source_elem.transition('zoom');
                    window.location.reload();
                    //window.location.href = data.redir;
                    popup_message(source_elem, "Moving Post", 'success', 1000);
                }
            },
            error: function (xhr, status, text) {
                error_message(elem, xhr, status, text);
                reset_drag(elem, source)
            }
        });

    ev.dataTransfer.clearData();
    dragged_over = '';

    reset_drag(elem, source)


}


function autocomplete_users(users_list) {
    // Add autocomplete to any text area element.
    var autocomplete = $('.autocomplete');

    // Map values in list to a list of dict [{key:'chosen', name:'displayed name'}....]
    var vals = $.map(users_list, function (value) {
        return {
            key: value,
            name: value
        };
    });

    function img_url(username) {
        let url = '/ajax/user/image/{0}/'.format(username);
        return "<li> <img class='ui circular image' style='display: inline' src={0} height='20' width='20' /> <b>{1}</b></li>".format(url, username);

    }

    // Autocomplete settings
    var AutocompleteSettings = {
        // Gets triggered at @
        at: "@",
        data: vals,
        displayTpl: img_url('${key}'),
        insertTpl: '@${key}',
        delay: 40
    };

    autocomplete.atwho(AutocompleteSettings);


}


function mark_spam(post_id, elem) {

    $.ajax('/ajax/report/spam/' + post_id + "/",
        {
            type: 'GET',
            dataType: 'json',
            ContentType: 'application/json',
            data: {},
            success: function (data) {
                //alert(elem.html());
                if (data.status === 'error') {
                    popup_message(elem.parent().parent(), data.msg, data.status);

                } else {
                    popup_message(elem.parent().parent(), data.msg, data.status);
                    $('#' + post_id).removeClass('open').addClass('spam');
                }

            },
            error: function (xhr, status, text) {
                error_message($('#' + post_id), xhr, status, text)
            }

        });


}


function moderate(elem, url) {

    event.preventDefault();
    //var elem = $(this);
    $('#modpanel').remove();

    // Could be a user or post uid
    var data_uid = elem.attr('data-value');

    var container = $('.moderate-insert[data-value="' + data_uid + '"]');
    var mod_url = url + data_uid + '/';
    //alert("FOOO");
    var page = $('<div id="modpanel"></div>').load(mod_url);
    container.after(page);
    //alert(container.html());
}

function similar_posts(elem) {
    var uid = elem.attr('post_uid');
    // Construct the similar posts link.
    var feed_url = '/similar/posts/' + uid + '/';
    var dimm_elem = $('#dim-similar');
    dimm_elem.dimmer('show');

    $.ajax(feed_url,
        {
            type: 'GET',
            dataType: 'json',
            ContentType: 'application/json',
            data: {
                'uid': uid
            },
            success: function (data) {
                if (data.status === 'error') {
                    popup_message(elem, data.msg, data.status);
                    dimm_elem.dimmer('hide');
                } else {
                    // Populate the feed.
                    dimm_elem.dimmer('hide');
                    elem.html(data.html);

                }
            },
            error: function (xhr, status, text) {
                error_message(elem, xhr, status, text)
            }
        })
}

function change_digest(elem, item) {
    var active = $('#digest-active');
    var icon_container = $('#digest-icon');
    var icon_str = item.data('icon');
    // Subscription url
    var digest_url = '/ajax/digest/';

    $.ajax(digest_url,
        {
            type: 'POST',
            dataType: 'json',
            ContentType: 'application/json',
            data: {
                'pref': value
            },
            success: function (data) {
                if (data.status === 'error') {
                    popup_message(elem, data.msg, data.status);
                } else {
                    // Replace current item with the select one.
                    active.text(item.text());
                    icon_container.removeClass();
                    icon_container.addClass(icon_str);
                }
            },
            error: function (xhr, status, text) {
                error_message(elem, xhr, status, text)
            }
        })
}

function change_subs(elem, value, $item) {
    var post_id = elem.attr("data-uid");
    // Currently selected item
    var active = $('#sub-active');
    // Subscription url
    var subs_url = '/ajax/subscribe/';

    $.ajax(subs_url,
        {
            type: 'POST',
            dataType: 'json',
            ContentType: 'application/json',
            data: {
                'root_uid': post_id,
                'sub_type': value
            },
            success: function (data) {
                if (data.status === 'error') {
                    popup_message(elem, data.msg, data.status);
                } else {
                    // Replace current item with the select one.
                    active.text($item.text());
                }

            },
            error: function (xhr, status, text) {
                error_message(elem, xhr, status, text)
            }
        })
}

$(document).ready(function () {


    $('.mark-spam.item').click(function (event) {
        var post_id = $(this).closest('.post').attr('id');

        //alert($(this).closest('.post').attr('id'));
        mark_spam(post_id, $(this));
    });


    function allowDrop(ev) {
        ev.preventDefault();
    }

    $('.post.comment').mousedown(function () {
        $(this).css('cursor', 'grabbing')
    });

    $('#similar-feed').each(function () {
        var elem = $(this);
        similar_posts(elem);
    });

    $('.ui.dropdown').dropdown();

    $('.editable').click(function (event) {
        if (event.metaKey || event.ctrlKey) {
            inplace_post_edit($(this))
        }
    }).dblclick(function (event) {
        inplace_post_edit($(this))
    });

    $('.inplace-click').click(function (event) {
        event.preventDefault();
        var uid = $(this).data('value');
        var elem = $('.editable[data-value="' + uid + '"]');
        inplace_post_edit(elem);
    });

    $(this).on('keyup', '#wmd-input', function (event) {
        var text = $(this).val();

        var preview = $('#preview');
        highlight(text, preview);

    });

    $(this).on('click', '#wmd-button-bar', function (event) {
        setTimeout(function () {
            var text = $('#wmd-input').val();
            var preview = $('#preview');
            highlight(text, preview);
        }, 10);

    });

    $(this).on('click', '#wmd-button-bar.answer', function (event) {
        setTimeout(function () {
            var text = $('textarea.answer').val();
            var preview = $('#answer-preview');
            highlight(text, preview);
        }, 10);

    });


    $(this).on('keyup', 'textarea.answer', function (event) {
        var text = $(this).val();
        var preview = $('#answer-preview');
        highlight(text, preview);
    });


    $(this).on('click', 'textarea.answer', function (event) {
        setTimeout(function () {
            var text = $('textarea.answer').val();
            var preview = $('#answer-preview');
            highlight(text, preview);
        }, 10);

    });


    $(this).keyup(function (event) {
        if (event.keyCode === 27) {
            $('.inplace').each(function () {
                event.preventDefault();
                var uid = $(this).data("value");
                cancel_inplace(uid);
                $('#new-comment').remove();
                $('#add-answer').html('');
            });
        }

    });


    $('#digest').dropdown({
        action: 'hide',
        onChange: function (value, text, $item) {
            var elem = $(this);
            change_digest(elem, $item);
        }
    });

    $('#subscribe')
        .dropdown({
            action: 'hide',
            onChange: function (value, text, $item) {
                var elem = $(this);
                change_subs(elem, value, $item);
            }
        });

    $(".add-comment").click(function (event) {
        event.preventDefault();

        var create_url = '/inplace/form/';
        var parent_uid = $(this).data('value');

        var container = $('.comment-insert[data-post="' + parent_uid + '"]');
        var post_actions = $('.hide-on-comment[data-value="' + parent_uid + '"]');

        cancel_inplace();

        // Check for existing comment.
        var comment = $("#new-comment");

        if (comment.length) {
            // Remove comment if exists.
            comment.remove();
        } else {
            // Create a new comment.
            comment = $('<div id="new-comment"></div>')
        }
        container.after(comment);
        comment.html('');
        //alert(container.length);

        $.ajax(create_url,
            {
                type: 'GET',
                dataType: 'json',
                ContentType: 'application/json',
                data: {
                    'uid': parent_uid,
                    'add_comment': 1,
                },
                success: function (data) {
                    if (data.status === 'error') {

                        popup_message(container, data.msg, data.status);

                    } else {
                        post_actions.hide();
                        comment.css({'padding-top': '5px', 'padding-bottom': '5px'});
                        comment.html(data.inplace_form);
                        //container.transition('slide down', 250);
                        comment.find('#wmd-input').focus();
                    }
                },
                error: function (xhr, status, text) {
                    error_message($(this), xhr, status, text)
                }
            })
    });

    remove_trigger();

    $(".moderate-post").click(function (event) {
        event.preventDefault();
        var elem = $(this);
        moderate(elem, '/moderate/');
    });

    $(".moderate-user").click(function (event) {
        event.preventDefault();
        var elem = $(this);
        moderate(elem, '/accounts/moderate/');

    });

    $('.vote').each(function (event) {
        var elem = $(this);
        var data_state = elem.attr('data-state');

        // Set the on class if the vote is selected.
        if (data_state === "1") {
            elem.addClass("on")
        }

        // Actions taken on vote click.
        $($(this)).click(function () {
            var elem = $(this);
            var post_uid = elem.attr('data-value');
            var data_type = elem.attr('data-type');
            apply_vote(elem, post_uid, data_type);
        });
    });

    $("#form-errors .error").each(function () {

        var elem = $(this);
        // Get errored out field id and label
        var field_id = elem.attr('data-value');
        var field_label = elem.attr('label');
        // Get the error message
        var message = elem.attr("message");
        // Select field in the form using it's id
        var field = $(field_id);
        // Add an 'error' to '.ui.field' to turn it red.
        field.closest(".field").addClass("error");
        // Insert the error message
        field.before('<div class="ui small red message"> {1}</div>'.f(field_label, message))
    });

    $(this).on('click', '#cancel', function () {
        $('#new-comment').remove();
        cancel_inplace();
    });

    $('.ui.sticky').sticky();
    //$('#chat-drop')
    $('pre').addClass('language-bash');
    $('code').addClass('language-bash');
    Prism.highlightAll();


    $('.tag-field').dropdown({

        allowAdditions: true,
        // Get form field to add to
        onChange: function (value, text, $selectedItem) {
            // Get form field to add to

            var tagid = $(this).children('select').attr('field_id');
            var tag_field = $('#{0}'.f(tagid));
            // Add selected tag to field
            tag_field.val(value);

        }
    });

    $('.tag-field > input.search').keydown(function (event) {
        // Prevent submitting form when adding tag by pressing ENTER.
        if (event.keyCode === 13) {
            event.preventDefault();
        }
        // Set value with SPACE bar
        if (event.keyCode === 32) {
            event.preventDefault();
            //alert(alert( $(this).parent('.tag-field').children('select').attr('class')));
            $(this).parent('.tag-field').children('select').dropdown('set selected', $(this).val().trim());
            $(this).val('');
        }

    });

    // $('.watched-tag-field').dropdown({
    //     allowAdditions: true,
    //     onChange: function (value, text, $selectedItem) {
    //         var tagid = $("#tag-menu").attr('field_id');
    //         var tag_field = $('#{0}'.f(tagid));
    //         // Add selected tag to field
    //         //alert(value);
    //         // Add selected tag to field
    //         //alert(value);
    //         tag_field.val(value);
    //     }
    // });
    //
    // $('.watched-tag-field >input.search').keydown(function (event) {
    //     // Prevent submitting form when adding tag by pressing ENTER.
    //     if (event.keyCode === 13) {
    //         event.preventDefault();
    //     }
    //     // Set value with SPACE bar
    //     if (event.keyCode === 32) {
    //         event.preventDefault();
    //         $("#watched-tags").dropdown('set selected', $(this).val().trim());
    //         $(this).val('')
    //
    //     }
    //
    // });

    $('#show-answer').click(function () {
        $('.hidden-answer').show()
    });

    $('.vote').popup({
        on: 'hover'
    });

    $('.trigger-highlight').each(function (event) {
        var elem = $(this);
        let container = $('#search-results');
        let query = container.data('query');
        let stop_words = container.data('stop');

        if (container.html() === undefined || container.html() === null) {
        } else {
            let stop_list = stop_words.split(',');
            highlight_search(query, elem, stop_list)
        }

    });

    var converter = new Markdown.getSanitizingConverter();
    var editor = new Markdown.Editor(converter);
    editor.run();
})
;
