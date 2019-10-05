function check_job() {


    //var job_uid = $(this).data('value');
    $('.check_back').each(function () {
    var job_uid = $(this).data('value');
    var state = $(this).data('state');
    if (job_uid === null || job_uid === undefined) {
        return
    }

    //alert('goo');
    $.ajax('/ajax/check/job/' + job_uid + '/',{
            type: 'GET',
            dataType: 'json',
            data:{'state': state},
            ContentType: 'application/json',
            success: function (data) {
                //alert($('.job-container-'+  job_uid).html());
                if (data.state_changed){
                    $('.job-container-' + job_uid).html(data.html);
                    var job_item = $('.job-item-' + job_uid);
                    job_item.transition('flash');
                }

            },
            error: function () {
            }
    })
    });


}

$(document).ready(function () {


     $('select')
        .dropdown()
    ;

//    $(".items > .item").click(function (event) {
//        var obj = $(this).find("a:first");
//        if (typeof obj !== 'undefined') {
//            window.location = obj.attr("href");
//       }
//    });
    setInterval(check_job, 20000);

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
