Example settings for various deployment options.


## Changing site domain

Inside of `conf/run/site_settings.py` edit the following settings:

    SITE_DOMAIN = "< production domain >"
    ALLOWED_HOSTS = [SITE_DOMAIN, ... ]


Inside of `conf/run/site_nginx.conf`, add the site domain to the server block :
    
    
    # Default server configuration.
    server {
        
        # DOMAIN NAME ADDED
        server_name < domain name >
        
        listen 80;
     }


## Change credentials 


### Email credentials
Inside of `conf/run/site_secrets.py` change the following email settings :
    
    DEFAULT_FROM_EMAIL = ""
    EMAIL_BACKEND = 'biostar.emailer.backend.SSLEmailBackend'
    EMAIL_HOST = 'email-smtp.us-west-2.amazonaws.com'
    EMAIL_HOST_USER = "< host user > "
    EMAIL_HOST_PASSWORD = "< host password > "
    EMAIL_PORT = 465

### Adding Oauth credentials 

Inside of `conf/run/site_secrets.py` edit the `SOCIAL_CLIENTS` list:


    SOCIAL_CLIENTS = [
    
         ("Name",
          "Client ID",
          "Client Secret"),
          
         ("GitHub",
          "id",
          "secret"),
      ]

**Note** : The server needs to be migrate to see the changes 

## Adding RECAPTCHA keys

Inside of `conf/run/site_secrets.py` edit the two settings:
 
    RECAPTCHA_PUBLIC_KEY = ""
    RECAPTCHA_PRIVATE_KEY = ""

 