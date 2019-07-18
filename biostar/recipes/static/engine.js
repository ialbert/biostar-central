


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
