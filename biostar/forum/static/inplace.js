function hidden_selector(elem) {

    // Hide or show the title, voting buttons, and content
    // when working with inplace forms.

    var title = ".post >  .title";
    var voting = ".post > .body > .voting";
    var content = ".post > .body > .content > .inplace";

    // Concatenate the final selector
    var selector = "{0},{1},{2}".format(title, voting, content);
    if (elem) {
        return elem.closest(selector)
    } else return $(selector)

}

function create_comment() {

    // Get the fields to submit
    var elem = $("#new-content");
    var parent = elem.closest('.post');
    var uid = parent.data("value");

    var content = $('#wmd-input');
    var cap_response = captcha();


    $.ajax('/ajax/comment/create/',
        {
            type: 'POST',
            dataType: 'json',
            ContentType: 'html',
            traditional: true,
            data: {
                'content': content.val(),
                'parent': uid,
                'recaptcha_response': cap_response,
            },
            success: function (data) {
                if (data.status === 'error') {
                    popup_message(parent, data.msg, data.status, 3000);
                    popup_message(elem, data.msg, data.status, 3000);
                } else {
                    // Redirect to the new post view.
                    window.location = data.redirect;
                    window.location.reload();
                }
            },
            error: function (xhr, status, text) {
                error_message(elem, xhr, status, text)
            }
        })
}

function cancel_inplace() {

    // Find the inplace element
    var inplace = $('#new-content');

    //if (inplace.length){}
    // Remove inplace item
    inplace.remove();

    // Find the hidden items
    var hidden = hidden_selector();

    // Show hidden items.
    hidden.each(function () {
        var is_hidden = $(this).is(':hidden');
        //
        if (is_hidden) {
            //alert($(this).attr('class'));
            hidden.show();
        }

    });


}


function prepare_inplace(content) {

    var post = content.closest(".post");

    // Get the current inplace items set to be replaced.
    var container = content.closest(".post > .body > .content > .inplace");

    // Any previously open inplace form gets closed.
    cancel_inplace(post);

    return container

}

function activate_inplace(current, form, hide, elem) {
    // TODO: need to rename variables, hard to understand.
    current.hide();
    current.html(form);
    current.show(300).find('textarea').focus();

    //var preview = $('#preview');
    current.find('pre').addClass('language-python');
    current.find('code').addClass('language-bash');

    if (hide) {
        //alert(hide);
        var to_hide = hidden_selector(elem);
        to_hide.hide();
    }


}


function inplace_form(elem, add_comment) {

    add_comment = add_comment || false;
    // Hiding the parent comment is dependent on whether a comment is being added
    var hide_containter = !add_comment;
    var post = elem.closest(".post");
    var uid = post.data("value");
    var url = '/inplace/form/';

    // Prepare the container
    var container = prepare_inplace(elem);
    var current = $('<div id="new-content"></div>');
    container.after(current);

    // Prepare the request data
    var input = {'uid': uid};
    if (add_comment) {
        input['add_comment'] = 1;
    } else {
        container.dimmer('show');
    }

    $.ajax(url,
        {
            type: 'GET',
            dataType: 'json',
            ContentType: 'application/json',
            data: input,
            success: function (data) {

                if (data.status === 'error') {
                    popup_message(post, data.msg, data.status, 3000);
                    return
                }
                activate_inplace(current, data.inplace_form, hide_containter, elem);
                container.dimmer('hide');

            },
            error: function (xhr, status, text) {
                error_message(elem, xhr, status, text)
            }
        })
}


function update_post(post, data) {

    // Current post content and title to replace
    // with returned values.
    var post_content = post.find('.editable');
    var post_title = post.find('.title');
    var post_tags = post.find('.tags');


    // Replace current post info with edited data
    post_content.html(data.html).show().focus();
    // Highlight text in content
    activate_prism(post_content);

    post_title.html(data.title).show();
    post_tags.html(data.tag_html).show();

    cancel_inplace(post);
}


function edit_post(post) {

    var uid = post.data('value');

    var edit_url = '/ajax/edit/{0}/'.format(uid);

    // Get the element being
    var form = $('#new-edit');

    // Form elements to submit.
    var title = form.find('#title');
    var content = form.find('#wmd-input');
    var type = form.find('#type').dropdown('get value');
    var tags = form.find('#tag-menu').dropdown('get value');

    title = title.val() || '';
    if (!($.isNumeric(type))) {
        type = -1
    }

    var cap_response = captcha();

    $.ajax(edit_url,
        {
            type: 'POST',
            dataType: 'json',
            ContentType: 'application/json',
            traditional: true,
            data: {
                'content': content.val(),
                'title': title,
                'type': type,
                'tag_val': tags,
                'recaptcha_response': cap_response
            },
            success: function (data) {
                if (data.status === 'error') {
                    popup_message(form, data.msg, data.status, 3000);
                } else {
                    // Update post with latest
                    update_post(post, data);
                }
            },
            error: function (xhr, status, text) {
                error_message(form, xhr, status, text)
            }
        })
}


$(document).on(function () {

    // Initialize pagedown

    $('.ui.dropdown').dropdown();

    $(this).on('click', '.editable', function (event) {
        if (event.metaKey || event.ctrlKey) {
            inplace_form($(this))
        }
    }).dblclick(function (event) {
        inplace_form($(this))
    });

    $(this).on('click', '.edit-button', function (event) {
        event.preventDefault();
        inplace_form($(this));
    });

    $(this).on('click', '.post .modal .exit.button', function () {
        var modal = $(this).closest('.modal');
        var post = $(this).closest('.post');
        modal.modal('hide');
        cancel_inplace(post);
    });

    $(this).on('click', '.post .modal .stay.button', function () {
        var modal = $(this).closest('.modal');
        modal.modal('hide');
    });

    $(this).on('click', '#inplace .cancel', function () {
        var post = $(this).closest('.post');
        cancel_inplace(post);
    });

    $(this).on('click', '#inplace .save', function () {
        event.preventDefault();
        var post = $(this).closest('.post');
        edit_post(post);
    });
    $(this).on('click', '#inplace .create', function () {
        event.preventDefault();
        create_comment();
    });

    $(".add-comment").click(function (event) {
        event.preventDefault();
        inplace_form($(this), true);

    });

    $(this).keyup(function (event) {
        if (event.keyCode === 27) {
            $('#inplace').each(function () {
                event.preventDefault();
                var post = $(this).closest('.post');
                cancel_inplace(post);
            });
        }
    });

    $('#wmd-input').keyup(function (event) {

        // Submit form with CTRL-ENTER
        if (event.ctrlKey && (event.keyCode === 13 || event.keyCode === 10)) {
            var save = $('#inplace').find('.save, .create');
            save.click();
        }
    });

});