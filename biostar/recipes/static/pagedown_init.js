var DjangoPagedown = DjangoPagedown | {};

DjangoPagedown = (function() {
  var converter = Markdown.getSanitizingConverter();
  var editors = {};
  var elements;

  Markdown.Extra.init(converter, {
    extensions: "all"
  });

  var createEditor = function(element) {
    var input = element.getElementsByClassName("wmd-input")[0];
    var id = input.id.substr(9);
    if (!editors.hasOwnProperty(id)) {
      var editor = new Markdown.Editor(converter, id, {});

      // Handle image upload
      if (element.classList.contains("image-upload-enabled")) {
        var upload = element.getElementsByClassName("pagedown-image-upload")[0];
        var url = upload.getElementsByClassName("url-input")[0];
        var file = upload.getElementsByClassName("file-input")[0];
        var cancel = upload.getElementsByClassName("deletelink")[0];
        var submit = upload.getElementsByClassName("submit-input")[0];

        var close = function(value, callback) {
          upload.classList.remove("show");
          url.value = "";
          file.value = "";
          callback(value);
        };

        editor.hooks.set("insertImageDialog", function(callback) {
          upload.classList.add("show");

          cancel.addEventListener(
            "click",
            function(event) {
              close(null, callback);
              event.preventDefault();
            },
            { once: true }
          );

          submit.addEventListener(
            "click",
            function() {
              // Regular URL
              if (url.value.length > 0) {
                close(url.value, callback);
              }
              // File upload
              else if (file.files.length > 0) {
                var data = new FormData();
                var xhr = new XMLHttpRequest();
                data.append("image", file.files[0]);
                xhr.open("POST", "/pagedown/image-upload/", true);
                xhr.addEventListener(
                  "load",
                  function() {
                    if (xhr.status !== 200) {
                      alert(xhr.statusText);
                    } else {
                      var response = JSON.parse(xhr.response);
                      if (response.success) {
                        close(response.url, callback);
                      } else {
                        if (response.error) {
                          var error = "";
                          for (var key in response.error) {
                            if (response.error.hasOwnProperty(key)) {
                              error += key + ":" + response.error[key];
                            }
                          }
                          alert(error);
                        }
                        close(null, callback);
                      }
                    }
                  },
                  {
                    once: true
                  }
                );
                xhr.send(data);
              } else {
                // Nothing
                close(null, callback);
              }
              event.preventDefault();
            },
            { once: true }
          );

          return true;
        });
      }

      editor.run();
      editors[id] = editor;
    }
  };

  var destroyEditor = function(element) {
    if (editors.hasOwnProperty(element.id)) {
      delete editors[el.id];
      return true;
    }
    return false;
  };

  var init = function() {
    elements = document.getElementsByClassName("wmd-wrapper");
    for (var i = 0; i < elements.length; ++i) {
      createEditor(elements[i]);
    }
  };

  return {
    init: function() {
      return init();
    },
    createEditor: function(element) {
      return createEditor(element);
    },
    destroyEditor: function(element) {
      return destroyEditor(element);
    }
  };
})();

window.onload = DjangoPagedown.init;
