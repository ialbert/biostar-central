
function create_comment(parent) {

    var create_url = '/ajax/comment/create/';
    // Get the fields to submit
    var form_elem = $('#post-form');
    var content = $('#wmd-input');

    var cap_response = captcha();


    $.ajax(create_url,
        {
            type: 'POST',
            dataType: 'json',
            ContentType: 'html',
            traditional: true,
            data: {
                'content': content.val(),
                'parent': '{{ post.uid }}',
                'recaptcha_response': cap_response,
            },
            success: function (data) {
                if (data.status === 'error') {
                    popup_message($('.error-msg'), data.msg, data.status, 3000);
                    popup_message(form_elem, data.msg, data.status, 3000);
                } else {
                    // Redirect user to the new post view.
                    window.location = data.redirect;
                    window.location.reload();
                }
            },
            error: function (xhr, status, text) {
                error_message(form_elem, xhr, status, text)
            }
        })
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

    $('#new-comment').remove();
    //Delete the form
    inplace_content.remove();
    //inplace_title.html("");
    // Hide the container
    // Show original content
    content.show();
    //Show any blocked element
    hidden.show();

}


function inplace_form(elem, add) {

    add = add || false;
    var post = elem.closest(".post");
    var url = '/inplace/form/';
    var uid = post.data('value');

    var container = post.find(".content div ");
    var hide = post.find(".title, .voting, .actions, .content div ");
    var current = $('<div id="new-content"></div>');

     // Any inplace forms already open get closed.
    cancel_inplace();
    container.dimmer('show');
    container.after(current);

    var input = {'uid': uid};

    if (add) {
        input.push({'add_comment': 1});
    }

    $.ajax(url,
        {
            type: 'GET',
            dataType: 'json',
            ContentType: 'application/json',
            data:input,
            success: function (data) {

                if (data.status === 'error') {
                    popup_message(elem, data.msg, data.status, 3000);
                    return
                }
                if (!add){
                    hide.hide();
                }
                current.html(data.inplace_form);
                current.show().find('textarea').focus();

                var preview = $('#preview');
                preview.find('pre').addClass('language-python');
                preview.find('code').addClass('language-bash');

                container.dimmer('hide');

            },
            error: function (xhr, status, text) {
                error_message(elem, xhr, status, text)
            }
        })
}

function edit_post(post, uid, elem) {

    var edit_url = '/ajax/edit/' + uid + '/';

    // Get form element
    var form_elem = elem.closest('form');

    // Inplace form container
    var form_container = $('#new-edit');
    // Hidden elements
    var hidden = $('.hide-on-edit');

    // Post title inside of the form
    var title = $('#title');
    var content = $('#wmd-input');
    var post_type = $('#type').dropdown('get value');
    var tag_val = $('#tag-menu').dropdown('get value');

    // Current post content and title to replace
    // with returned values.
    var post_content = $('.editable[data-value="' + uid + '"]');
    var post_title = $('.post-title[data-value="' + uid + '"]');
    var post_tags = $('.post-tags[data-value="' + uid + '"]');

    title = title.val() || '';
    if (!($.isNumeric(post_type))) {
        post_type = -1
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
                'type': post_type,
                'tag_val': tag_val,
                'recaptcha_response': cap_response
            },
            success: function (data) {
                if (data.status === 'error') {
                    popup_message($('.error-msg'), data.msg, data.status, 3000);
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
                    //highlight()
                    activate_prism(post_content);
                }
            },
            error: function (xhr, status, text) {
                error_message(form_elem, xhr, status, text)
            }
        })
}

function add_comment(parent, elem) {


    var create_url = '/inplace/form/';
    var parent_uid = elem.data('value');

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

}


$(document).on(function () {

    // Initialize pagedown
    init_pagedown();

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

    $(this).on('click', '#inplace .modal .exit.button', function () {
        var modal = $(this).closest('.modal');
        modal.modal('hide');
        cancel_inplace();
    });

    $(this).on('click', '#inplace .modal .stay.button', function () {
        var modal = $(this).closest('.modal');
        modal.modal('hide');
    });

    $(this).on('click', '#inplace .cancel', function () {
        cancel_inplace();
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

    $(this).on('click', ".add-comment", function (event) {
        event.preventDefault();
        inplace_form($(this), true)
    });

    $(this).keyup(function (event) {
        if (event.keyCode === 27) {
            $('.inplace').each(function () {
                event.preventDefault();
                var uid = $(this).data("value");
                cancel_inplace(uid);
            });
        }
    });

    $('#wmd-input').keyup(function (event) {

        // Submit form with CTRL-ENTER
        if (event.ctrlKey && (event.keyCode === 13 || event.keyCode === 10)) {
            var save = $('#inplace').find('.save,.create');
            save.click();
        }
    });


});