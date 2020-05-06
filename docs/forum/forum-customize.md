# Documentation for Customizing the Biostars Forum


## Setting up remote postgres database for synchronization 

**Note: You need SSH access to the remote server.**

### 1. Login to remote server.


     ssh user@host      # Log into remote server

### 2. Get location of `postgresql.conf`.
 
 
    find / -name "postgresql.conf" # Find the postgresql.conf
    
The location is expected to look like : `/etc/postgresql/10/main/postgresql.conf`

### 3. Edit `postgresql.conf`. 

Open postgresql.conf file and add the following line to the end.

     listen_addresses = '*'

This allows for postgres to listen to all IP address on the remote server. Without this, it is restricted to 127.0.0.1 


### 4. Find and edit `pg_hba.conf`

Repeat step 2 to find `pg_hba.conf` file:

     find / -name "pg_hba.conf"

It is expected to be in the same directory as the previously edited  `postgresql.conf`.

The location is expected to look like : `/etc/postgresql/10/main/pg_hba.conf`

Add the following line to the end of `pg_hba.conf` file after replacing the `IP`:

    # IP = remote IP from which connection is allowed
    
    host    all         all     < IP >/32   md5  

`md5` stands for the type of authentication type used when making a TCP connection. In this case, it is an md5 hash. 

For more on PostgresSQL authentication: https://www.postgresql.org/docs/9.6/auth-pg-hba-conf.html


### 5. Restart PostgreSQL server to apply the changes.


    # Restart method 1
    service postgresql restart
    
    # Restart method 2
    sudo pg_ctlcluster 10 main start 
