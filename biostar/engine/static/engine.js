$(document).ready(function () {

    $('select.dropdown')
        .dropdown()
    ;

    $(".item").click(function (event) {
        var obj = $(this).find("a:first");
        window.location = obj.attr("href");
    });

});
