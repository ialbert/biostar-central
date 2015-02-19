$(document).ready(function () {

    // Set focus on search field.
    var search = $('#search');

    if (search.length) {
        size = search.val().length;
        search.focus();
        search[0].setSelectionRange(size, size);
    }

    // Draw pietyis.
    $(".line").peity("line")

});
