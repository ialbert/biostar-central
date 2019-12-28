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
    elem.toggleClass("on")

    // The vote handler.
    vote_url = "/ajax/vote/"

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


function edit_post(uid) {

    var edit_url = '/ajax/edit/' + uid + '/';
    // Rendered form element
    var form_elem = $('#post-form');
    // Inplace form container
    var form_container = $('#new-edit');
    // Hidden elements
    var hidden = $('.hide-on-edit');

    // Post title inside of the form
    var title = $('#title');
    var content = $('#wmd-input');
    var post_type = $('#type').dropdown('get value');
    var tag_val = $('#tags').dropdown('get value');

    // Current post content and title to replace
    // with returned values.
    var post_content = $('.editable[data-value="' + uid + '"]');
    var post_title = $('.post-title[data-value="' + uid + '"]');
    var post_tags = $('.post-tags[data-value="' + uid + '"]');

    //alert(content.val());
    // Title is null and type are null
    // meaning current post is not top level.
    if (title.val() == null) {
        title.val('')
    }
    if (!($.isNumeric(post_type))) {
        post_type = -1
    }

    $.ajax(edit_url,
        {
            type: 'POST',
            dataType: 'json',
            ContentType: 'application/json',
            traditional: true,
            data: {
                'content': content.val(),
                'title': title.val(),
                'type': post_type,
                'tag_val': tag_val,
            },
            success: function (data) {
                if (data.status === 'error') {
                    popup_message(form_elem, data.msg, data.status, 3000);
                } else {

                    // Clear and hide inplace form
                    form_elem.html('');
                    form_container.hide();

                    // Show hidden items
                    hidden.show();

                    // Replace current post info with edited data
                    post_content.html(data.html).show().focus();
                    post_title.html(data.title).show();
                    post_tags.html(data.tag_html).show();

                    // Highlight the markdown in content
                    post_content.find('pre').addClass('language-bash');
                    post_content.find('code').addClass('language-bash');
                    Prism.highlightAll();
                }
            },
            error: function (xhr, status, text) {
                error_message(form_elem, xhr, status, text)
            }
        })
}


function cancel_create() {

    var form_container = $('#insert-form');
    form_container.html('');
    $('.new-post').removeClass('active');
    $('.dim-on-create').removeClass('fit-create');

}

function cancel_inplace() {

    var inplace_content = $('#new-edit');
    //var inplace_title = $('inplace-title[data-value="'+ uid +'"]');

    //var title = $('.editable-title[data-value="'+ uid +'"]');
    var content = $('.editable');
    //var content = $('#wmd-input-' + uid);
    var hidden = $('.hide-on-edit');
    $('.hide-on-comment').show();

    $('.dim-on-create').dimmer('hide').removeClass('fit-create');
    ;
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
        if ($.inArray(value, stop_list) === -1){
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
    cancel_create();
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
                } else {
                    elem.hide();
                    hidden.hide();
                    dim_elem.dimmer('hide');
                    editing.html(data.inplace_form);
                    editing.show().find('textarea').focus();
                    var preview = $('#preview');
                    preview.find('pre').addClass(' language-python');
                    preview.find('code').addClass('  language-bash language-python');
                    //preview.find('code').removeClass('lang-python');
                    //Prism.highlightAll();

                }
            },
            error: function (xhr, status, text) {
                error_message(elem, xhr, status, text)
            }
        })
}


function search(query, elem, search_url) {
    var res = $('#results');

    var container = $('#contain-search');
    var redir = container.data('redir');

    container.addClass('loading search');
    res.width(container.width());
    res.addClass('ui search message');
    res.html('Searching ...');

    $.ajax(search_url, {
        type: 'GET',
        dataType: 'json',
        ContentType: 'application/json',
        data: {
            'query': query,
            'redir': redir,
        },

        success: function (data) {
            if (data.status === 'error') {
                popup_message($(this), data.msg, data.status);
            } else {
                // Success
                //alert(data.html);
                res.removeClass('ui message');
                res.html('');

                if (data.html === null || data.html === undefined) {
                    window.location.href = data.redir
                }
                res.html(data.html);
                container.removeClass('loading search');
            }

        },
        error: function (xhr, status, text) {
            error_message(elem, xhr, status, text)
        }
    });
}


function highlight(text, preview) {

    var con = markdownit();
    var res = con.render(text);
    //var preview = $('#preview');
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
    let elem = $('.droptarget[data-value="' + elem_id + '"]');

    children = '{0}'.format(elem.data('thread'));

    ev.dataTransfer.setData("text", elem_id);
    ev.target.style.opacity = "0.4";

    //alert(elem_id);

}

function drag_leave(ev, pid) {
    let elem = $('.droptarget[data-value="' + pid + '"]');
    elem.css('border', 'none');

    //elem.css('padding', '0');
    elem.css('opacity', '1');
    //elem.css('backgroundColor', 'none');
    dragged_over = '';
    //elem.removeClass('foo');
    //ev.stopPropagation();
    //elem.removeClass('dragging-over');

}


function drag_over(ev, pid) {
    //ev.preventDefault();
    //var source = ev.dataTransfer.getData("text");
    let elem = $('.droptarget[data-value="' + pid + '"]');
    let children_list = children.split(",");
//elem.addClass('foo');
    dragged_over = pid;

    if (jQuery.inArray(pid, children_list) !== -1) {
        //elem.css('backgroundColor', '#ffb5a8');
        elem.css('border', '#ffb5a8 dashed');

        //dragged_over = '';
        //alert(children_list)
    } else {

        //alert(pid);
        elem.css('border', '#c2ffc2 dashed 5px');
        //alert("FOOO");
        //elem.css('bac')
        //elem.css('padding', '1px');
        //elem.css('backgroundColor', '#c2ffc2')
    }
    //elem.addClass('dragging-over');
//ev.stopPropagation();
    //ev.stopPropagation();


}


function drop(ev, elem_id) {

    ev.preventDefault();

    var source = ev.dataTransfer.getData("text");
    let source_elem = $('#' + source);
    source_elem.css('opacity', '1');
    let elem = $('.droptarget[data-value="' + elem_id + '"]');

    let children_list = children.split(",");

    //let children_list = children.split(",");

    //alert(children_list);
    //alert(elem_id);

    if (jQuery.inArray(elem_id, children_list) === -1 && dragged_over === elem_id) {

        //alert(dragged_over);
        $.ajax('/drag/and/drop/',
            {
                type: 'POST',
                dataType: 'json',
                ContentType: 'application/json',
                data: {
                    'uid': source,
                    'parent': elem_id,

                },
                success: function (data) {
                    if (data.status === 'error') {
                        popup_message(elem, data.msg, data.status);
                        elem.css('border', 'none');
                        source_elem.css('border', 'none');
                        elem.css('opacity', '1');
                        elem.removeClass('foo');
//elem.removeClass('dragging-over');
                    } else {
                        //alert(data.redir)
                        //window.location.href = data.redir;
                        //window.location.reload();

                        //elem.removeClass('foo');
                        source_elem.transition('zoom');
                        //elem.style.backgroundColor = '';
                        //ev.target.style.border = "red";
                        //elem.removeChild = document.getElementById('indent-' + source)
                        //elem.append(document.getElementById('indent-' + source));

                        //window.location.href = data.redir;
                        window.location.reload();
                        popup_message(source_elem, "Moving Post", 'success', 1000);
                        //source_elem.parent().transition('pulse');
                        //source_elem.parent().transition('pulse');
                    }
                },
                error: function (xhr, status, text) {
                    error_message(elem, xhr, status, text);
                    elem.css('border', 'none');
                    source_elem.css('border', 'none');
                    elem.css('opacity', '1');
                }
            });

        //elem.css('border', '');
        //elem.css('paddingRight', '#0');
        //elem.css('backgroundColor', 'white');

        ev.dataTransfer.clearData();
        dragged_over = '';
        //ev.stopPropagation();
    } else if (dragged_over === elem_id) {

        popup_message(source_elem, "Can not be dropped here.", 'error', 1000);

    }
    elem.css('border', 'none');
    source_elem.css('border', 'none');
    elem.css('opacity', '1');


}

function chat_list(user_uid) {
    // Show list of active chats a user has
    var chat_list_url = '/ajax/chat/list/';
    var container = $('#contain-chat');

    $.ajax(chat_list_url,
        {
            type: 'GET',
            dataType: 'json',
            ContentType: 'application/json',
            data: {
                'uid': user_uid
            },
            success: function (data) {
                if (data.status === 'error') {
                    popup_message(container, data.msg, data.status);
                } else {
                    //var contain = "<div class='ui message'>{0}, most recent: {1}</div>".f(data.nposts, data.most_recent_url);
                    container.html(data.html)
                }
            },
            error: function (xhr, status, text) {
                error_message(container, xhr, status, text)
            }
        })

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

$(document).ready(function () {


    $('#chat').click(function () {
        var user_uid = $(this).data('value');
        chat_list(user_uid);
    });

    $(this).click(function (event) {
        var res = $('#results');
        if (typeof res.html() === 'undefined') {
            return
        }
        if (res.html().length > 0) {
            res.html('');
            res.removeClass('ui message');
        }
    });

    function allowDrop(ev) {
        ev.preventDefault();
    }

    $('.post.comment').mousedown(function () {
        $(this).css('cursor', 'grabbing')
    })

    $('#similar-feed').each(function () {
        var elem = $(this);
        //elem.html('asshole');
        // Fetch feed from url url
        var uid = elem.attr('post_uid');
        var feed_url = '/similar/posts/' + uid + '/';


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
                    } else {
                        // Populate the feed.
                        elem.html(data.html)
                    }
                },
                error: function (xhr, status, text) {
                    error_message(elem, xhr, status, text)
                }
            })
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
        var uid = $(this).data('value');
        var elem = $('.editable[data-value="' + uid + '"]');
        inplace_post_edit(elem);
    });

    $(this).on('keyup', '#wmd-input', function (event) {
        var text = $(this).val();

        var preview = $('#preview');
        highlight(text, preview);

    });

    $(this).on('keydown', '#wmd-input-id_content', function (event) {
        var text = $(this).val();
        var preview1 = $('#wmd-preview-wmd-input-id_content');

        highlight(text, preview1);

        var preview = $('#wmd-preview-id_content');
        highlight(text, preview);
    });


    $(this).on('click', '#wmd-button-bar-wmd-input-id_content', function (event) {
        setTimeout(function () {
            var text = $('#wmd-input-id_content').val();
            var preview = $('#wmd-preview-wmd-input-id_content');
            highlight(text, preview);
        }, 10);

    });
    $(this).on('click', '#wmd-button-bar-id_content', function (event) {
        setTimeout(function () {
            var text = $('#wmd-input-id_content').val();
            var preview = $('#wmd-preview-id_content');
            highlight(text, preview);
        }, 10);

    });


    $(this).on('click', '#wmd-button-bar', function (event) {
        setTimeout(function () {
            var text = $('#wmd-input').val();
            var preview = $('#preview');
            highlight(text, preview);
        }, 10);

    });

    $(this).keyup(function (event) {
        if (event.keyCode === 27) {
            $('.inplace').each(function () {
                event.preventDefault();
                var uid = $(this).data("value");
                cancel_inplace(uid);
                cancel_create();
                $('#new-comment').remove();
                $('#add-answer').html('');
            });
        }

    });


    $('#digest').dropdown({
        action: 'hide',
        onChange: function (value, text, $item) {
            var elem = $(this);

            console.log(elem.id);

            // Get the root id
            // Currently selected item
            var active = $('#digest-active');
            var icon_container = $('#digest-icon');
            var icon_str = $item.data('icon');
            // Subscription url
            var digest_url = '/ajax/digest/';
            //alert('ffff')
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
                            active.text($item.text());
                            icon_container.removeClass();
                            icon_container.addClass(icon_str);
                        }
                    },
                    error: function (xhr, status, text) {
                        error_message(elem, xhr, status, text)
                    }
                })
        }
    });

    $('#subscribe')
        .dropdown({
            action: 'hide',
            onChange: function (value, text, $item) {
                var elem = $(this);

                console.log(elem.id);

                // Get the root id
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
        });

    $(".add-comment").click(function (event) {
        event.preventDefault();

        var create_url = '/inplace/form/';
        var parent_uid = $(this).data('value');

        //var container = $("#comment-insert-" + parent_uid);
        var container = $("#comment-insert-" + parent_uid);
        //var url = "/new/comment/" + post_uid + "/";
        var post_actions = $('.hide-on-comment[data-value="' + parent_uid + '"]');
        cancel_inplace();
        $('#insert-form').html('');
        $('.new-post').removeClass('active');

        $('#add-answer').html('');

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
                    'parent': parent_uid,
                    'comment': 1,
                    'top': 0,
                    'rows': 6,
                },
                success: function (data) {
                    if (data.status === 'error') {
                        popup_message($('#error'), data.msg, data.status);
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

    $('#search').keyup(function (event) {
        var query = $(this).val();
        var search_url = $(this).attr('url');
        // Only preform searches when pressing ENTER
        if (event.keyCode === 13) {
            search(query, $(this), search_url);
        }
    });

    $('.new-post').click(function () {
        var create_url = '/inplace/form/';
        var form_container = $('#insert-form');
        cancel_inplace();
        $('#new-comment').remove();
        $('#add-answer').html('');
        $('.dim-on-create').dimmer('hide');

        $.ajax(create_url,
            {
                type: 'GET',
                dataType: 'json',
                ContentType: 'application/json',
                data: {
                    'rows': 13,
                },
                success: function (data) {
                    if (data.status === 'error') {
                        popup_message($('#error'), data.msg, data.status);

                    } else {
                        $('#menu-header > .item').each(function () {
                            $(this).removeClass('active');
                        });
                        $('.new-post').addClass('active');

                        form_container.html(data.inplace_form);
                        form_container.show();
                        form_container.find('#title').focus();
                        $('.dim-on-create').dimmer('show').addClass('fit-create');

                    }

                },
                error: function (xhr, status, text) {
                    error_message($(this), xhr, status, text)
                }
            })
    });

    $(".moderate-post").click(function (event) {
        event.preventDefault();
        var elem = $(this);

        $('#modpanel').remove();

        // Could be a user or post uid
        var data_uid = elem.attr('data-value');

        var container = $("#moderate-insert-" + data_uid);
        var mod_url = '/moderate/' + data_uid + '/';

        var page = $('<div id="modpanel"></div>').load(mod_url);
        container.after(page)

    });

    $('a').click(function () {
        $(this).transition('pulse', 295);
    });

    $(this).on('click', '.show-preview', function () {
        var contain = $('#html-preview');
        var preview = $('#preview');
        contain.transition('slide down', 400);
        preview.find('pre').addClass('language-bash');
        preview.find('code').addClass('language-bash');
        //Prism.highlightAll();
        //Prism.highlightAll(preview.find('code'));

    });

    $(".moderate-user").click(function (event) {
        //alert("FOO");
        event.preventDefault();
        var elem = $(this);

        $('#modpanel').remove();

        // Could be a user or post uid
        var data_uid = elem.attr('data-value');

        var container = $("#moderate-insert-" + data_uid);
        var mod_url = '/accounts/moderate/' + data_uid + '/';
//alert(container.html());

        var page = $('<div id="modpanel"></div>').load(escape(mod_url));
        container.after(page);

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
        cancel_create();
        cancel_inplace();
        $('#add-answer').html('');
    });
    $('.ui.sticky').sticky()
    ;

    $('.display-answer').click(function () {
        var create_url = '/inplace/form/';
        var form_container = $('#add-answer');
        var parent_uid = form_container.data('value');
        //var form_is_visible = $('#insert-form.visible').val() != null;

        $('#new-comment').remove();
        cancel_inplace();
        cancel_create();
        // Avoid hitting the server when the form is already visible.

        $.ajax(create_url,
            {
                type: 'GET',
                dataType: 'json',
                ContentType: 'application/json',
                data: {
                    'rows': 15,
                    'parent': parent_uid,
                    'top': 0,
                },
                success: function (data) {
                    if (data.status === 'error') {
                        popup_message($('#error'), data.msg, data.status);
                    } else {
                        form_container.html(data.inplace_form);
                        form_container.show();
                        form_container.find('#wmd-input').focus();
                    }

                },
                error: function (xhr, status, text) {
                    error_message($(this), xhr, status, text)
                }
            })
    });

    //$('#chat-drop')
    $('pre').addClass('language-bash');
    $('code').addClass('language-bash');
    Prism.highlightAll();
    var preview = $('#wmd-preview');

    $('#wmd-input').keyup(function () {

        preview.find('code').addClass('language-bash');
        preview.find('pre').addClass('language-bash');

    }).keydown(function () {
        preview.find('code').addClass('language-bash');
        preview.find('pre').addClass('language-bash');
    });

    $('#wmd-button-bar').click(function () {
        setTimeout(function () {
            preview.find('code').addClass('language-bash');
            preview.find('pre').addClass('language-bash');
        }, 10)
    });

    $('.tag-field').dropdown({

        allowAdditions: true,
        // Get form field to add to
        onChange: function (value, text, $selectedItem) {
            // Get form field to add to
            var tagid = $("#tag-menu").attr('field_id');
            var tag_field = $('#{0}'.f(tagid));
            // Add selected tag to field
            //alert(value);
            tag_field.val(value);
        }
    });

    $('.tag-field >input.search').keydown(function (event) {

        //event.stopPropagation();
        // Prevent submitting form when adding tag by pressing ENTER.
        if (event.keyCode === 13) {
            event.preventDefault();
            //alert($('#my_tags_id').val());
            //event.stopPropagation();
        }
        // Set value with SPACE bar
        if (event.keyCode === 32) {
            event.preventDefault();
            //event.stopPropagation();
            $("#tag-menu").dropdown('set selected', $(this).val().trim());
            $(this).val('');
            //alert($('#my_tags_id').val());
        }

    });

    $('.watched-tag-field').dropdown({
        allowAdditions: true,
        onChange: function (value, text, $selectedItem) {
            // Get form field to add to
            var tagid = $("#watched-tags").attr('field_id');
            var tag_field = $('#{0}'.f(tagid));
            // Add selected tag to field
            //alert(value);
            tag_field.val(value);
        }
    });

    $('.watched-tag-field >input.search').keydown(function (event) {
        //event.stopPropagation();
        // Prevent submitting form when adding tag by pressing ENTER.
        if (event.keyCode === 13) {
            event.preventDefault();
            //event.stopPropagation();
        }
        // Set value with SPACE bar
        if (event.keyCode === 32) {
            event.preventDefault();
            //event.stopPropagation();
            $("#watched-tags").dropdown('set selected', $(this).val().trim());
            $(this).val('')

        }

    });

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

        if (container.html() === undefined || container.html() === null){
        }else{
            let stop_list = stop_words.split(',');
            highlight_search(query, elem, stop_list)
        }

    });

    var autocomplete_elem = $('#wmd-input-id_content');

    if (autocomplete_elem.data('users') !== undefined || autocomplete_elem.data('users') != null) {
        var users = autocomplete_elem.data('users').split(',');
        autocomplete_users(users)
    }


})
;
