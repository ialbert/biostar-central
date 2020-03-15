
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
                    // Redirect user to the new post view.
                    window.location = data.redirect;
                    window.location.reload();
                }
            },
            error: function (xhr, status, text) {
                error_message(elem, xhr, status, text)
            }
        })
}

function cancel_inplace(post) {

    var inplace = $('#new-content');
    //var inplace_title = $('inplace-title[data-value="'+ uid +'"]');
    //var title = $('.editable-title[data-value="'+ uid +'"]');
    //var content = $('.editable');
    //var content = $('#wmd-input-' + uid);
    var hidden = post.find('.editable, .title, .voting, .actions, .content div ');
    //$('.hide-on-comment').show();

    //$('#new-comment').remove();
    //Delete the form
    inplace.remove();
    //inplace_title.html("");
    // Hide the container
    // Show original content
    //content.show();
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
    cancel_inplace(post);
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


function update_post(post, data){

    // Current post content and title to replace
    // with returned values.
    var post_content = post.find('.editable');
    var post_title = post.find('.title');
    var post_tags = post.find('.tags');

    // Replace current post info with edited data
    post_content.html(data.html).show().focus();
    post_title.html(data.title).show();
    post_tags.html(data.tag_html).show();

    cancel_inplace(post);
}



function edit_post(post) {

    var uid = post.data('value');

    var edit_url = '/ajax/edit/{0}/'.format(uid);

    // Get the element being
    var form = $('#new-edit');

    // Post title inside of the form
    var title = form.find('#title');
    var content = form.find('#wmd-input');
    var post_type = form.find('#type').dropdown('get value');
    var tag_val = form.find('#tag-menu').dropdown('get value');

    // Current post content and title to replace
    // with returned values.
    var post_content = post.find('.editable');

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
                    // Update post with lateets
                    update_post(post, data);

                    // Highlight text in content
                    activate_prism(post_content);
                }
            },
            error: function (xhr, status, text) {
                error_message(form_elem, xhr, status, text)
            }
        })
}


$(document).on(function () {

    // Initialize pagedown
    init_pagedown();
    tags_dropdown();
    activate_prism();

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

    $(this).on('click', '.modal .exit.button', function () {
        var modal = $(this).closest('.modal');
        var post = $(this).closest('.post');
        modal.modal('hide');
        cancel_inplace(post);
    });

    $(this).on('click', '.modal .stay.button', function () {
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

    $(this).on('click', ".add-comment", function (event) {
        event.preventDefault();
        inplace_form($(this), true)
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
            var save = $('#inplace').find('.save,.create');
            save.click();
        }
    });

});