$(document).ready(function () {

    // Set focus on search field.
    var search = $('#search');
    var size = search.val().length;
    search.focus();
    search[0].setSelectionRange(size, size);

});
