
$(document).ready(function () {

     $('select')
        .dropdown()
    ;

    $(".items > .item").click(function (event) {
        var obj = $(this).find("a:first");
        if (typeof obj !== 'undefined') {
            window.location = obj.attr("href");
        }
    });

});
