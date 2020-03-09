//function

function project_id() {
    return $("#project_id").val();

}

function change_state(job, data) {
    // Get the current job image, stdout and stderr
    var image = job.find('#img:first');
    var stdout = $('#stdout pre:first');
    var stderr = $('#stderr');

    // Only replace state and activate effects when the state has actually changed.
    if (data.state_changed) {
        job.find('.state:first').html(data.html);
        job.transition('pulse');

    }
    // Update the image, stdout, and stderr.
    image.replaceWith(data.img_tmpl);
    stdout.text(data.stdout);
    stderr.text(data.stderr);

}

function trigger_running(job, data) {

    var loader = $('#stdout .loader');
    var uid = job.data('value');
    var link = job.find(".link");

    if (data.is_running) {
        // Anchor link to the log when job is running
        link.attr("href", '/job/view/{0}/#log'.format(uid));
        // Add the "Running" loader under log
        loader.html('<div id="log" class="ui log message">' +
            '<span class="ui active small inline loader"></span>' +
            '<span>Running</span>' +
            '</div>');
    } else {
        loader.html("");
        link.attr("href", '/job/view/{0}/'.format(uid))
    }
}


function check_jobs() {

    // Look at each job with a 'check_back' tag
    $('.check_back').each(function () {

        var job = $(this).closest('.job');

        // Get the uid and state
        var uid = job.data('value');

        // Get the current job state.
        var state = $(this).data("state");

        // Bail out when a job uid is not provided.
        if (uid === null || uid === undefined) {
            return
        }

        $.ajax('/ajax/check/job/{0}/'.format(uid), {
            type: 'GET',
            dataType: 'json',
            data: {'state': state},
            ContentType: 'application/json',
            success: function (data) {

                // Change the job state
                change_state(job, data);

                // Trigger events when job is running
                trigger_running(job, data);

                // Redirect to the file list once the job is done
                // and the current page is a job view.
                if (data.redir && $("#view").length) {
                    window.location.replace(data.redir + "#flist");
                    window.location.reload()
                }
            },
            error: function (xhr, status, text) {
                error_message($(this), xhr, status, text)
            }
        })
    });

}


function copy_object(uid, clipboard, container) {

    $.ajax('/copy/object/',
        {
            type: 'POST',
            dataType: 'json',
            data: {'uid': uid, 'clipboard': clipboard},

            success: function (data) {
                if (data.status === 'success') {
                    // Get the item container
                    container.transition({
                        animation: 'pulse', onComplete: function () {
                            container.addClass('copied item')
                        }
                    });
                    update_clipboard();
                    return
                }
                popup_message(container, data.msg, data.status, 4000)

            },
            error: function (xhr, status, text) {
                error_message(container, xhr, status, text)
            }


        }
    )
}

function copy_file(path, elem) {

    $.ajax('/file/copy/', {
            type: 'POST',
            dataType: 'json',
            data: {'path': path, 'uid': project_id()},
            success: function (data) {
                if (data.status === 'success') {
                    popup_message(elem, data.msg, data.status, 500);
                    update_clipboard();
                    return
                }
                popup_message(elem, data.msg, data.status, 2000)

            },
            error: function (xhr, status, text) {
                error_message(elem, xhr, status, text)
            }
        }
    )
}


function paste(data){
    // Paste a given item into the clipboard
    var elem = $("#clipboard");
    $.ajax('/paste/', {
            type: 'POST',
            dataType: 'json',
            data:data,
            success: function (data) {
                if (data.status === "success") {

                     window.location.href = data.redirect;
                     popup_message(elem, data.msg, data.status, 500);
                }
                //popup_message(elem, data.msg, data.status, 2000)
            },
            error: function (xhr, status, text) {
                error_message(elem, xhr, status, text)
            }
        }
    )
}

function clear(){
    var elem = $("#clipboard");
    $.ajax('/clear/', {
            type: 'POST',
            dataType: 'json',
            success: function (data) {
                if (data.status === "success") {
                     popup_message(elem, data.msg, data.status, 500);
                     elem.hide();
                     elem.html(data.html);
                     elem.fadeIn("slow", function () {});
                     $('.copied.item').each(function () {
                         $(this).removeClass(" copied")
                     })
                    return
                }
                popup_message(elem, data.msg, data.status, 4000)
            },
            error: function (xhr, status, text) {
                error_message(elem, xhr, status, text)
            }
        }
    )
}

function update_clipboard(){
    var elem = $("#clipboard");
    $.ajax('/clipboard/', {
            type: 'POST',
            dataType: 'json',
            data:{"id": project_id()},
            success: function (data) {
                if (data.html) {
                     elem.hide();
                     elem.html(data.html);
                     elem.fadeIn("slow", function () {});
                }
                //popup_message(elem, data.msg, data.status, 2000)
            },
            error: function (xhr, status, text) {
                error_message(elem, xhr, status, text)
            }
        }
    )
}

function toggle_delete(elem, otype) {

    var uid = elem.data("value");
    // Fetch the element to decrement/increment  once toggle is complete.
    var count = $("#{0}_count".format(otype));
    var url = '/toggle/delete/';
    var data = {'uid': uid, 'type': otype};

    $.ajax(url, {
            type: 'POST',
            dataType: 'json',
            data: data,
            success: function (data) {
                if (data.status === 'success') {
                    // Give message
                    elem.transition('zoom');
                    // Update counts on tabular menu
                    count.html(data.counts);
                    return
                }
                popover_message(elem.before(), data.msg, data.status, 2000)
            },
            error: function (xhr, status, text) {
                error_message(elem.before(), xhr, status, text)
            }
        }
    )


}

function change_access(access, user_id, project_uid, elem) {
    let container = $('.container-' + user_id);
    $.ajax('/manage/access/', {
            type: 'POST',
            dataType: 'json',
            data: {
                'user_id': user_id,
                'access': access,
                'project_uid': project_uid
            },

            success: function (data) {
                if (data.status === 'success') {
                    container.transition('bounce').transition({
                        animation: 'pulse', onComplete: function () {
                            if (data.no_access) {
                                container.removeClass('inputcolor')
                            } else {
                                container.addClass('inputcolor')
                            }
                        }
                    });
                    return
                }
                popup_message(container, data.msg, data.status, 2000)

            },
            error: function (xhr, status, text) {
                error_message(elem, xhr, status, text)
            }
        }
    )

}

$(document).ready(function () {

    $('.access-dropdown').dropdown();
    $(this).on('click', '.access-dropdown .menu .item', function () {
        let user = $(this).parent().data('user');
        let project = $(this).parent().data('project');
        let value = $(this).data("value");
        change_access(value, user, project)

    });

    $('.ui.dropdown').dropdown();

    $('select').dropdown();

    setInterval(check_jobs, 5000);

    $('#recipe-search').keyup(function () {
        var query = $(this).val();
        $.ajax("/search/", {
            type: 'GET',
            dataType: 'html',
            data: {'q': query},
            success: function (data) {
                $('#search-results').html(data)
            },
            error: function () {
            }
        });

    });


    $('.listing').popup({
        on: 'hover',
    });

    $(".moderate-user").click(function (event) {
        event.preventDefault();
        var elem = $(this);

        $('#modpanel').remove();

        // Could be a user or post uid
        var data_uid = elem.attr('data-value');

        var container = $("#moderate-insert-" + data_uid);
        var mod_url = '/accounts/moderate/' + data_uid + '/';

        var page = $('<div id="modpanel"></div>').load(mod_url);
        container.after(page)
    });

    $(this).on('click', '.job .delete', function (event) {
        event.preventDefault();
        var elem = $(this).closest('.job');
        toggle_delete(elem, 'job')
    });

    $(this).on('click', '.data .delete', function (event) {
        event.preventDefault();
        var elem = $(this).closest('.data');
        toggle_delete(elem, 'data')
    });


    $('.checkbox').checkbox();

    $(this).on('click', '.data .copy.button', function () {
        let data = $(this).closest('.data');
        let uid = data.data('value');
        copy_object(uid, "data", data);
    });

    $(this).on('click', '.job .copy.button', function () {
        let job = $(this).closest('.job');
        let uid = job.data('value');
        copy_object(uid, "job", job);
    });

    $(this).on('click', '.recipe .copy.button', function () {
        let recipe = $(this).closest('.recipe');
        let uid = recipe.data("value");
        copy_object(uid, "recipe", recipe);
    });

    $(this).on('click', '.file .copy', function () {
        let file = $(this).closest('.file');
        let path = file.data("value");
        copy_file(path, file);
    });

    $('#clipboard').each(function () {
        update_clipboard();
    });

    $(this).on('click', '#clipboard .paste', function () {
        var data = {"id": project_id()};
        paste(data);
    });

    $(this).on('click', '#clipboard .clone', function () {
        var data = {"id": project_id(), "target":"clone"};
        paste(data);
    });

    $(this).on('click', '#clipboard .clear', function () {
        clear();
    });

    $('pre').addClass('language-bash');
    $('code').addClass('language-bash').css('padding', '0');
    Prism.highlightAll();

});
