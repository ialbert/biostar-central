





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

        event.preventDefault();
        var elem = $(this);
        var data_uid = elem.attr('data-uid');

        $.ajax("/data/copy", {
                type: 'GET',
                dataType: 'json',
                data: {data_uid: data_uid},
                success: function (data) {
                    alert(data.message);
                        },
                });

    });

});
