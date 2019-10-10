
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


function check_job() {

    // Look at each job with a 'check_back' tag
    $('.check_back').each(function () {
        // Get the uid and state
        var job_uid = $(this).data('value');
        var state = $(this).data('state');
        // Bail out when a job uid is not provided.
        if (job_uid === null || job_uid === undefined) {
            return
        }

        // Ajax request checking state change and replace appropriate element.
        $.ajax('/ajax/check/job/' + job_uid + '/', {
            type: 'GET',
            dataType: 'json',
            data: {'state': state},
            ContentType: 'application/json',
            success: function (data) {
                // Only replace and activate effects when the state has actually changed.
                if (data.state_changed) {
                    $('.job-container-' + job_uid).html(data.html);
                    var job_item = $('.job-item-' + job_uid);
                    job_item.transition('pulse');
                }
            },
            error: function (xhr, status, text) {
            error_message($(this), xhr, status, text)
        }
        })
    });

}

$(document).ready(function () {

    $('.ui.dropdown').dropdown();
     $('select')
        .dropdown()
    ;
     $('.ui.sticky')
  .sticky()
;


     $('#preview').click(function (event) {



     });
      $('#json_preview').click(function (event) {
         let recipe_uid = $(this).data('value');
         alert(recipe_uid)


     });

//    $(".items > .item").click(function (event) {
//        var obj = $(this).find("a:first");
//        if (typeof obj !== 'undefined') {
//            window.location = obj.attr("href");
//       }
//    });

    // Check and update 'Running' and 'Spooled' jobs every 20 seconds.
    setInterval(check_job, 5000 );

    $(".copy-data").click(function (event) {

        var elem = $(this);
        var data_uid = elem.attr('data-uid');
        var copy_url = elem.attr('copy-url');

        $.ajax(copy_url, {
                type: 'GET',
                dataType: 'json',
                ContentType: 'application/json',
                data: {data_uid: data_uid},
                success: function (data) {
                pop_over($("#copy-message-"+ data_uid), data.msg, data.status );
                },
                error: function () {
                }

                });

    });

    $('#recipe-search').keyup(function(){

        var query = $(this).val();

        $.ajax("/search/", {
            type: 'GET',
            dataType: 'html',

            data: {'q': query},

            success: function (data) {
            $('#search-results').html(data);
            },
            error: function () {
            }
            });

    });

});
