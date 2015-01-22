PRAGMA foreign_keys=OFF;
BEGIN TRANSACTION;
CREATE TABLE "auth_permission" (
    "id" integer NOT NULL PRIMARY KEY,
    "name" varchar(50) NOT NULL,
    "content_type_id" integer NOT NULL,
    "codename" varchar(100) NOT NULL,
    UNIQUE ("content_type_id", "codename")
);
INSERT INTO "auth_permission" VALUES(1,'Can add permission',1,'add_permission');
INSERT INTO "auth_permission" VALUES(2,'Can change permission',1,'change_permission');
INSERT INTO "auth_permission" VALUES(3,'Can delete permission',1,'delete_permission');
INSERT INTO "auth_permission" VALUES(4,'Can add group',2,'add_group');
INSERT INTO "auth_permission" VALUES(5,'Can change group',2,'change_group');
INSERT INTO "auth_permission" VALUES(6,'Can delete group',2,'delete_group');
INSERT INTO "auth_permission" VALUES(7,'Can add content type',3,'add_contenttype');
INSERT INTO "auth_permission" VALUES(8,'Can change content type',3,'change_contenttype');
INSERT INTO "auth_permission" VALUES(9,'Can delete content type',3,'delete_contenttype');
INSERT INTO "auth_permission" VALUES(10,'Can add site',4,'add_site');
INSERT INTO "auth_permission" VALUES(11,'Can change site',4,'change_site');
INSERT INTO "auth_permission" VALUES(12,'Can delete site',4,'delete_site');
INSERT INTO "auth_permission" VALUES(13,'Can add message body',5,'add_messagebody');
INSERT INTO "auth_permission" VALUES(14,'Can change message body',5,'change_messagebody');
INSERT INTO "auth_permission" VALUES(15,'Can delete message body',5,'delete_messagebody');
INSERT INTO "auth_permission" VALUES(16,'Can add message',6,'add_message');
INSERT INTO "auth_permission" VALUES(17,'Can change message',6,'change_message');
INSERT INTO "auth_permission" VALUES(18,'Can delete message',6,'delete_message');
INSERT INTO "auth_permission" VALUES(19,'Can add log entry',7,'add_logentry');
INSERT INTO "auth_permission" VALUES(20,'Can change log entry',7,'change_logentry');
INSERT INTO "auth_permission" VALUES(21,'Can delete log entry',7,'delete_logentry');
INSERT INTO "auth_permission" VALUES(22,'Can add flat page',8,'add_flatpage');
INSERT INTO "auth_permission" VALUES(23,'Can change flat page',8,'change_flatpage');
INSERT INTO "auth_permission" VALUES(24,'Can delete flat page',8,'delete_flatpage');
INSERT INTO "auth_permission" VALUES(25,'Can add session',9,'add_session');
INSERT INTO "auth_permission" VALUES(26,'Can change session',9,'change_session');
INSERT INTO "auth_permission" VALUES(27,'Can delete session',9,'delete_session');
INSERT INTO "auth_permission" VALUES(28,'Can add email address',10,'add_emailaddress');
INSERT INTO "auth_permission" VALUES(29,'Can change email address',10,'change_emailaddress');
INSERT INTO "auth_permission" VALUES(30,'Can delete email address',10,'delete_emailaddress');
INSERT INTO "auth_permission" VALUES(31,'Can add email confirmation',11,'add_emailconfirmation');
INSERT INTO "auth_permission" VALUES(32,'Can change email confirmation',11,'change_emailconfirmation');
INSERT INTO "auth_permission" VALUES(33,'Can delete email confirmation',11,'delete_emailconfirmation');
INSERT INTO "auth_permission" VALUES(34,'Can add migration history',12,'add_migrationhistory');
INSERT INTO "auth_permission" VALUES(35,'Can change migration history',12,'change_migrationhistory');
INSERT INTO "auth_permission" VALUES(36,'Can delete migration history',12,'delete_migrationhistory');
INSERT INTO "auth_permission" VALUES(37,'Can add user',13,'add_user');
INSERT INTO "auth_permission" VALUES(38,'Can change user',13,'change_user');
INSERT INTO "auth_permission" VALUES(39,'Can delete user',13,'delete_user');
INSERT INTO "auth_permission" VALUES(40,'Can add email list',14,'add_emaillist');
INSERT INTO "auth_permission" VALUES(41,'Can change email list',14,'change_emaillist');
INSERT INTO "auth_permission" VALUES(42,'Can delete email list',14,'delete_emaillist');
INSERT INTO "auth_permission" VALUES(43,'Can add tag',15,'add_tag');
INSERT INTO "auth_permission" VALUES(44,'Can change tag',15,'change_tag');
INSERT INTO "auth_permission" VALUES(45,'Can delete tag',15,'delete_tag');
INSERT INTO "auth_permission" VALUES(46,'Can add profile',16,'add_profile');
INSERT INTO "auth_permission" VALUES(47,'Can change profile',16,'change_profile');
INSERT INTO "auth_permission" VALUES(48,'Can delete profile',16,'delete_profile');
INSERT INTO "auth_permission" VALUES(49,'Can add tag',17,'add_tag');
INSERT INTO "auth_permission" VALUES(50,'Can change tag',17,'change_tag');
INSERT INTO "auth_permission" VALUES(51,'Can delete tag',17,'delete_tag');
INSERT INTO "auth_permission" VALUES(52,'Can add post',18,'add_post');
INSERT INTO "auth_permission" VALUES(53,'Can change post',18,'change_post');
INSERT INTO "auth_permission" VALUES(54,'Can delete post',18,'delete_post');
INSERT INTO "auth_permission" VALUES(55,'Can add reply token',19,'add_replytoken');
INSERT INTO "auth_permission" VALUES(56,'Can change reply token',19,'change_replytoken');
INSERT INTO "auth_permission" VALUES(57,'Can delete reply token',19,'delete_replytoken');
INSERT INTO "auth_permission" VALUES(58,'Can add email sub',20,'add_emailsub');
INSERT INTO "auth_permission" VALUES(59,'Can change email sub',20,'change_emailsub');
INSERT INTO "auth_permission" VALUES(60,'Can delete email sub',20,'delete_emailsub');
INSERT INTO "auth_permission" VALUES(61,'Can add email entry',21,'add_emailentry');
INSERT INTO "auth_permission" VALUES(62,'Can change email entry',21,'change_emailentry');
INSERT INTO "auth_permission" VALUES(63,'Can delete email entry',21,'delete_emailentry');
INSERT INTO "auth_permission" VALUES(64,'Can add post view',22,'add_postview');
INSERT INTO "auth_permission" VALUES(65,'Can change post view',22,'change_postview');
INSERT INTO "auth_permission" VALUES(66,'Can delete post view',22,'delete_postview');
INSERT INTO "auth_permission" VALUES(67,'Can add vote',23,'add_vote');
INSERT INTO "auth_permission" VALUES(68,'Can change vote',23,'change_vote');
INSERT INTO "auth_permission" VALUES(69,'Can delete vote',23,'delete_vote');
INSERT INTO "auth_permission" VALUES(70,'Can add subscription',24,'add_subscription');
INSERT INTO "auth_permission" VALUES(71,'Can change subscription',24,'change_subscription');
INSERT INTO "auth_permission" VALUES(72,'Can delete subscription',24,'delete_subscription');
INSERT INTO "auth_permission" VALUES(73,'Can add badge',25,'add_badge');
INSERT INTO "auth_permission" VALUES(74,'Can change badge',25,'change_badge');
INSERT INTO "auth_permission" VALUES(75,'Can delete badge',25,'delete_badge');
INSERT INTO "auth_permission" VALUES(76,'Can add award',26,'add_award');
INSERT INTO "auth_permission" VALUES(77,'Can change award',26,'change_award');
INSERT INTO "auth_permission" VALUES(78,'Can delete award',26,'delete_award');
INSERT INTO "auth_permission" VALUES(79,'Can add blog',27,'add_blog');
INSERT INTO "auth_permission" VALUES(80,'Can change blog',27,'change_blog');
INSERT INTO "auth_permission" VALUES(81,'Can delete blog',27,'delete_blog');
INSERT INTO "auth_permission" VALUES(82,'Can add blog post',28,'add_blogpost');
INSERT INTO "auth_permission" VALUES(83,'Can change blog post',28,'change_blogpost');
INSERT INTO "auth_permission" VALUES(84,'Can delete blog post',28,'delete_blogpost');
INSERT INTO "auth_permission" VALUES(85,'Can add social app',29,'add_socialapp');
INSERT INTO "auth_permission" VALUES(86,'Can change social app',29,'change_socialapp');
INSERT INTO "auth_permission" VALUES(87,'Can delete social app',29,'delete_socialapp');
INSERT INTO "auth_permission" VALUES(88,'Can add social account',30,'add_socialaccount');
INSERT INTO "auth_permission" VALUES(89,'Can change social account',30,'change_socialaccount');
INSERT INTO "auth_permission" VALUES(90,'Can delete social account',30,'delete_socialaccount');
INSERT INTO "auth_permission" VALUES(91,'Can add social token',31,'add_socialtoken');
INSERT INTO "auth_permission" VALUES(92,'Can change social token',31,'change_socialtoken');
INSERT INTO "auth_permission" VALUES(93,'Can delete social token',31,'delete_socialtoken');
INSERT INTO "auth_permission" VALUES(94,'Can add task state',32,'add_taskmeta');
INSERT INTO "auth_permission" VALUES(95,'Can change task state',32,'change_taskmeta');
INSERT INTO "auth_permission" VALUES(96,'Can delete task state',32,'delete_taskmeta');
INSERT INTO "auth_permission" VALUES(97,'Can add saved group result',33,'add_tasksetmeta');
INSERT INTO "auth_permission" VALUES(98,'Can change saved group result',33,'change_tasksetmeta');
INSERT INTO "auth_permission" VALUES(99,'Can delete saved group result',33,'delete_tasksetmeta');
INSERT INTO "auth_permission" VALUES(100,'Can add interval',34,'add_intervalschedule');
INSERT INTO "auth_permission" VALUES(101,'Can change interval',34,'change_intervalschedule');
INSERT INTO "auth_permission" VALUES(102,'Can delete interval',34,'delete_intervalschedule');
INSERT INTO "auth_permission" VALUES(103,'Can add crontab',35,'add_crontabschedule');
INSERT INTO "auth_permission" VALUES(104,'Can change crontab',35,'change_crontabschedule');
INSERT INTO "auth_permission" VALUES(105,'Can delete crontab',35,'delete_crontabschedule');
INSERT INTO "auth_permission" VALUES(106,'Can add periodic tasks',36,'add_periodictasks');
INSERT INTO "auth_permission" VALUES(107,'Can change periodic tasks',36,'change_periodictasks');
INSERT INTO "auth_permission" VALUES(108,'Can delete periodic tasks',36,'delete_periodictasks');
INSERT INTO "auth_permission" VALUES(109,'Can add periodic task',37,'add_periodictask');
INSERT INTO "auth_permission" VALUES(110,'Can change periodic task',37,'change_periodictask');
INSERT INTO "auth_permission" VALUES(111,'Can delete periodic task',37,'delete_periodictask');
INSERT INTO "auth_permission" VALUES(112,'Can add worker',38,'add_workerstate');
INSERT INTO "auth_permission" VALUES(113,'Can change worker',38,'change_workerstate');
INSERT INTO "auth_permission" VALUES(114,'Can delete worker',38,'delete_workerstate');
INSERT INTO "auth_permission" VALUES(115,'Can add task',39,'add_taskstate');
INSERT INTO "auth_permission" VALUES(116,'Can change task',39,'change_taskstate');
INSERT INTO "auth_permission" VALUES(117,'Can delete task',39,'delete_taskstate');
INSERT INTO "auth_permission" VALUES(118,'Can add queue',40,'add_queue');
INSERT INTO "auth_permission" VALUES(119,'Can change queue',40,'change_queue');
INSERT INTO "auth_permission" VALUES(120,'Can delete queue',40,'delete_queue');
INSERT INTO "auth_permission" VALUES(121,'Can add message',41,'add_message');
INSERT INTO "auth_permission" VALUES(122,'Can change message',41,'change_message');
INSERT INTO "auth_permission" VALUES(123,'Can delete message',41,'delete_message');
CREATE TABLE "auth_group_permissions" (
    "id" integer NOT NULL PRIMARY KEY,
    "group_id" integer NOT NULL,
    "permission_id" integer NOT NULL REFERENCES "auth_permission" ("id"),
    UNIQUE ("group_id", "permission_id")
);
CREATE TABLE "auth_group" (
    "id" integer NOT NULL PRIMARY KEY,
    "name" varchar(80) NOT NULL UNIQUE
);
CREATE TABLE "django_content_type" (
    "id" integer NOT NULL PRIMARY KEY,
    "name" varchar(100) NOT NULL,
    "app_label" varchar(100) NOT NULL,
    "model" varchar(100) NOT NULL,
    UNIQUE ("app_label", "model")
);
INSERT INTO "django_content_type" VALUES(1,'permission','auth','permission');
INSERT INTO "django_content_type" VALUES(2,'group','auth','group');
INSERT INTO "django_content_type" VALUES(3,'content type','contenttypes','contenttype');
INSERT INTO "django_content_type" VALUES(4,'site','sites','site');
INSERT INTO "django_content_type" VALUES(5,'message body','messages','messagebody');
INSERT INTO "django_content_type" VALUES(6,'message','messages','message');
INSERT INTO "django_content_type" VALUES(7,'log entry','admin','logentry');
INSERT INTO "django_content_type" VALUES(8,'flat page','flatpages','flatpage');
INSERT INTO "django_content_type" VALUES(9,'session','sessions','session');
INSERT INTO "django_content_type" VALUES(10,'email address','account','emailaddress');
INSERT INTO "django_content_type" VALUES(11,'email confirmation','account','emailconfirmation');
INSERT INTO "django_content_type" VALUES(12,'migration history','south','migrationhistory');
INSERT INTO "django_content_type" VALUES(13,'user','users','user');
INSERT INTO "django_content_type" VALUES(14,'email list','users','emaillist');
INSERT INTO "django_content_type" VALUES(15,'tag','users','tag');
INSERT INTO "django_content_type" VALUES(16,'profile','users','profile');
INSERT INTO "django_content_type" VALUES(17,'tag','posts','tag');
INSERT INTO "django_content_type" VALUES(18,'post','posts','post');
INSERT INTO "django_content_type" VALUES(19,'reply token','posts','replytoken');
INSERT INTO "django_content_type" VALUES(20,'email sub','posts','emailsub');
INSERT INTO "django_content_type" VALUES(21,'email entry','posts','emailentry');
INSERT INTO "django_content_type" VALUES(22,'post view','posts','postview');
INSERT INTO "django_content_type" VALUES(23,'vote','posts','vote');
INSERT INTO "django_content_type" VALUES(24,'subscription','posts','subscription');
INSERT INTO "django_content_type" VALUES(25,'badge','badges','badge');
INSERT INTO "django_content_type" VALUES(26,'award','badges','award');
INSERT INTO "django_content_type" VALUES(27,'blog','planet','blog');
INSERT INTO "django_content_type" VALUES(28,'blog post','planet','blogpost');
INSERT INTO "django_content_type" VALUES(29,'social app','socialaccount','socialapp');
INSERT INTO "django_content_type" VALUES(30,'social account','socialaccount','socialaccount');
INSERT INTO "django_content_type" VALUES(31,'social token','socialaccount','socialtoken');
INSERT INTO "django_content_type" VALUES(32,'task state','djcelery','taskmeta');
INSERT INTO "django_content_type" VALUES(33,'saved group result','djcelery','tasksetmeta');
INSERT INTO "django_content_type" VALUES(34,'interval','djcelery','intervalschedule');
INSERT INTO "django_content_type" VALUES(35,'crontab','djcelery','crontabschedule');
INSERT INTO "django_content_type" VALUES(36,'periodic tasks','djcelery','periodictasks');
INSERT INTO "django_content_type" VALUES(37,'periodic task','djcelery','periodictask');
INSERT INTO "django_content_type" VALUES(38,'worker','djcelery','workerstate');
INSERT INTO "django_content_type" VALUES(39,'task','djcelery','taskstate');
INSERT INTO "django_content_type" VALUES(40,'queue','django','queue');
INSERT INTO "django_content_type" VALUES(41,'message','django','message');
CREATE TABLE "django_site" (
    "id" integer NOT NULL PRIMARY KEY,
    "domain" varchar(100) NOT NULL,
    "name" varchar(50) NOT NULL
);
INSERT INTO "django_site" VALUES(1,'localhost:8080','Site Name');
CREATE TABLE "messages_messagebody" (
    "id" integer NOT NULL PRIMARY KEY,
    "text" text NOT NULL,
    "author_id" integer NOT NULL,
    "subject" varchar(120) NOT NULL,
    "parent_msg_id" integer REFERENCES "messages_messagebody" ("id"),
    "sent_at" datetime NOT NULL
);
INSERT INTO "messages_messagebody" VALUES(1,'Hello <b>Biostar Community</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:11.531000');
INSERT INTO "messages_messagebody" VALUES(2,'Hello <b>Istvan Albert</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.464000');
INSERT INTO "messages_messagebody" VALUES(3,'Hello <b>Fabio</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.506000');
INSERT INTO "messages_messagebody" VALUES(4,'Hello <b>Jason</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.516000');
INSERT INTO "messages_messagebody" VALUES(5,'Hello <b>Zhenhai Zhang</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.528000');
INSERT INTO "messages_messagebody" VALUES(6,'Hello <b>Tom Koerber</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.539000');
INSERT INTO "messages_messagebody" VALUES(7,'Hello <b>Suk211</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.549000');
INSERT INTO "messages_messagebody" VALUES(8,'Hello <b>Lemon</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.562000');
INSERT INTO "messages_messagebody" VALUES(9,'Hello <b>Wubin Qu</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.572000');
INSERT INTO "messages_messagebody" VALUES(10,'Hello <b>Question Bot</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.598000');
INSERT INTO "messages_messagebody" VALUES(11,'Hello <b>Reka Albert</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.608000');
INSERT INTO "messages_messagebody" VALUES(12,'Hello <b>Yang Yang</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.619000');
INSERT INTO "messages_messagebody" VALUES(13,'Hello <b>Gue Su Chang</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.630000');
INSERT INTO "messages_messagebody" VALUES(14,'Hello <b>Zhaorong</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.640000');
INSERT INTO "messages_messagebody" VALUES(15,'Hello <b>Nickey</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.652000');
INSERT INTO "messages_messagebody" VALUES(16,'Hello <b>Renee</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.661000');
INSERT INTO "messages_messagebody" VALUES(17,'Hello <b>Yu</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.672000');
INSERT INTO "messages_messagebody" VALUES(18,'Hello <b>Gue Su</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.683000');
INSERT INTO "messages_messagebody" VALUES(19,'Hello <b>Mohammed Islaih</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.693000');
INSERT INTO "messages_messagebody" VALUES(20,'Hello <b>Alex Reynolds</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.704000');
INSERT INTO "messages_messagebody" VALUES(21,'Hello <b>User 4824</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.714000');
INSERT INTO "messages_messagebody" VALUES(22,'Hello <b>User 8226</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.724000');
INSERT INTO "messages_messagebody" VALUES(23,'Hello <b>Giovanni M Dall&#39;Olio</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.735000');
INSERT INTO "messages_messagebody" VALUES(24,'Hello <b>Etal</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.746000');
INSERT INTO "messages_messagebody" VALUES(25,'Hello <b>Fabio</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.756000');
INSERT INTO "messages_messagebody" VALUES(26,'Hello <b>Nicojo</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.766000');
INSERT INTO "messages_messagebody" VALUES(27,'Hello <b>Allen Yu</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.777000');
INSERT INTO "messages_messagebody" VALUES(28,'Hello <b>Curious George</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.788000');
INSERT INTO "messages_messagebody" VALUES(29,'Hello <b>User 5034</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.799000');
INSERT INTO "messages_messagebody" VALUES(30,'Hello <b>Pierre Lindenbaum</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.810000');
INSERT INTO "messages_messagebody" VALUES(31,'Hello <b>Marcos De Carvalho</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.823000');
INSERT INTO "messages_messagebody" VALUES(32,'Hello <b>Gustavo Costa</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.836000');
INSERT INTO "messages_messagebody" VALUES(33,'Hello <b>User 6996</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.848000');
INSERT INTO "messages_messagebody" VALUES(34,'Hello <b>Paul J. Davis</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.859000');
INSERT INTO "messages_messagebody" VALUES(35,'Hello <b>David Nusinow</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.869000');
INSERT INTO "messages_messagebody" VALUES(36,'Hello <b>Brentp</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.880000');
INSERT INTO "messages_messagebody" VALUES(37,'Hello <b>Razor</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.894000');
INSERT INTO "messages_messagebody" VALUES(38,'Hello <b>Simon Cockell</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.905000');
INSERT INTO "messages_messagebody" VALUES(39,'Hello <b>Konrad</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.918000');
INSERT INTO "messages_messagebody" VALUES(40,'Hello <b>Biorelated</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.928000');
INSERT INTO "messages_messagebody" VALUES(41,'Hello <b>Yann Abraham</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.939000');
INSERT INTO "messages_messagebody" VALUES(42,'Hello <b>Fernando Muñiz</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.953000');
INSERT INTO "messages_messagebody" VALUES(43,'Hello <b>User 1402</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.963000');
INSERT INTO "messages_messagebody" VALUES(44,'Hello <b>Andrew</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.974000');
INSERT INTO "messages_messagebody" VALUES(45,'Hello <b>Alex</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.985000');
INSERT INTO "messages_messagebody" VALUES(46,'Hello <b>Owen</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:12.995000');
INSERT INTO "messages_messagebody" VALUES(47,'Hello <b>Chris</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.006000');
INSERT INTO "messages_messagebody" VALUES(48,'Hello <b>Kelly O.</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.017000');
INSERT INTO "messages_messagebody" VALUES(49,'Hello <b>User 3704</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.028000');
INSERT INTO "messages_messagebody" VALUES(50,'Hello <b>Greg</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.039000');
INSERT INTO "messages_messagebody" VALUES(51,'Hello <b>Pedrobeltrao</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.055000');
INSERT INTO "messages_messagebody" VALUES(52,'Hello <b>Liam Thompson</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.065000');
INSERT INTO "messages_messagebody" VALUES(53,'Hello <b>Michael Barton</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.077000');
INSERT INTO "messages_messagebody" VALUES(54,'Hello <b>Schrodinger&#39;S Cat</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.090000');
INSERT INTO "messages_messagebody" VALUES(55,'Hello <b>Michael Dondrup</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.101000');
INSERT INTO "messages_messagebody" VALUES(56,'Hello <b>Brad Chapman</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.114000');
INSERT INTO "messages_messagebody" VALUES(57,'Hello <b>Paulati</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.127000');
INSERT INTO "messages_messagebody" VALUES(58,'Hello <b>Dave Bridges</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.139000');
INSERT INTO "messages_messagebody" VALUES(59,'Hello <b>Daniel Swan</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.150000');
INSERT INTO "messages_messagebody" VALUES(60,'Hello <b>Lukfor</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.165000');
INSERT INTO "messages_messagebody" VALUES(61,'Hello <b>Chris Fields</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.176000');
INSERT INTO "messages_messagebody" VALUES(62,'Hello <b>Darked89</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.189000');
INSERT INTO "messages_messagebody" VALUES(63,'Hello <b>Lmartinho</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.201000');
INSERT INTO "messages_messagebody" VALUES(64,'Hello <b>Manuel Corpas</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.212000');
INSERT INTO "messages_messagebody" VALUES(65,'Hello <b>Abhishek Tiwari</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.239000');
INSERT INTO "messages_messagebody" VALUES(66,'Hello <b>Neilfws</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.259000');
INSERT INTO "messages_messagebody" VALUES(67,'Hello <b>Piotr Byzia</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.273000');
INSERT INTO "messages_messagebody" VALUES(68,'Hello <b>Jeroen Van Goey</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.288000');
INSERT INTO "messages_messagebody" VALUES(69,'Hello <b>Vince</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.301000');
INSERT INTO "messages_messagebody" VALUES(70,'Hello <b>Tom Walsh</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.315000');
INSERT INTO "messages_messagebody" VALUES(71,'Hello <b>Egon Willighagen</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.339000');
INSERT INTO "messages_messagebody" VALUES(72,'Hello <b>Mndoci</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.354000');
INSERT INTO "messages_messagebody" VALUES(73,'Hello <b>Jeremy Leipzig</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.370000');
INSERT INTO "messages_messagebody" VALUES(74,'Hello <b>Mmarchin</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.385000');
INSERT INTO "messages_messagebody" VALUES(75,'Hello <b>Paolo</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.397000');
INSERT INTO "messages_messagebody" VALUES(76,'Hello <b>Andrew</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.409000');
INSERT INTO "messages_messagebody" VALUES(77,'Hello <b>Eleanor Howe</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.420000');
INSERT INTO "messages_messagebody" VALUES(78,'Hello <b>Jboveda</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.440000');
INSERT INTO "messages_messagebody" VALUES(79,'Hello <b>Tim</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.451000');
INSERT INTO "messages_messagebody" VALUES(80,'Hello <b>Daniel Jurczak</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.463000');
INSERT INTO "messages_messagebody" VALUES(81,'Hello <b>User 1073</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.473000');
INSERT INTO "messages_messagebody" VALUES(82,'Hello <b>Geoffjentry</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.489000');
INSERT INTO "messages_messagebody" VALUES(83,'Hello <b>User 3505</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.500000');
INSERT INTO "messages_messagebody" VALUES(84,'Hello <b>Anshu</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.511000');
INSERT INTO "messages_messagebody" VALUES(85,'Hello <b>Bryan Maloney</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.522000');
INSERT INTO "messages_messagebody" VALUES(86,'Hello <b>Andrew Su</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.535000');
INSERT INTO "messages_messagebody" VALUES(87,'Hello <b>Khader Shameer</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.546000');
INSERT INTO "messages_messagebody" VALUES(88,'Hello <b>Pierre</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.565000');
INSERT INTO "messages_messagebody" VALUES(89,'Hello <b>Yuri</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.580000');
INSERT INTO "messages_messagebody" VALUES(90,'Hello <b>Allpowerde</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.593000');
INSERT INTO "messages_messagebody" VALUES(91,'Hello <b>Mikael Huss</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.619000');
INSERT INTO "messages_messagebody" VALUES(92,'Hello <b>Perry</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.637000');
INSERT INTO "messages_messagebody" VALUES(93,'Hello <b>User 5106</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.655000');
INSERT INTO "messages_messagebody" VALUES(94,'Hello <b>Nir London</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.674000');
INSERT INTO "messages_messagebody" VALUES(95,'Hello <b>Bioinfo</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.688000');
INSERT INTO "messages_messagebody" VALUES(96,'Hello <b>Jc</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.700000');
INSERT INTO "messages_messagebody" VALUES(97,'Hello <b>Michael Hoffman</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.725000');
INSERT INTO "messages_messagebody" VALUES(98,'Hello <b>Nancy Parmalee</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.737000');
INSERT INTO "messages_messagebody" VALUES(99,'Hello <b>Avilella</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.761000');
INSERT INTO "messages_messagebody" VALUES(100,'Hello <b>Ryan</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-04-29 15:02:13.774000');
INSERT INTO "messages_messagebody" VALUES(101,'<a href="/u/2/">Istvan Albert</a> wrote on

<a href="/p/1/">Site Use Guidelines</a> : 
Here are a few guidelines:


The site&#39;s goal is to answer bioinformatics and systems biology related questions
Answer questions to gain reputation. 
Don&#39;t forget to vote for answers that you like! Registered users may vote on answers.
If you are the one asking the original question you may also sele
',2,'Site Use Guidelines',NULL,'2009-09-30 20:12:07.053000');
INSERT INTO "messages_messagebody" VALUES(102,'<a href="/u/2/">Istvan Albert</a> wrote on

<a href="/p/2/">How Do I Convert From Bed Format To Gff Format?</a> : 
I have a file in GFF format and I need to convert it to BED format. What do I do?

',2,'How Do I Convert From Bed Format To Gff Format?',NULL,'2009-09-30 20:55:00.133000');
INSERT INTO "messages_messagebody" VALUES(103,'<a href="/u/2/">Istvan Albert</a> wrote on

<a href="/p/2/#3">A: How Do I Convert From Bed Format To Gff Format?</a> : 
Both formats are tab delimited text files used to represent DNA features in genomes. The order of columns between the two are different, there are also columns that correspond to attributes missing from one or the other format. Nonetheless the most important difference between the two is the coordin
',2,'A: How Do I Convert From Bed Format To Gff Format?',NULL,'2009-09-30 20:56:18.353000');
INSERT INTO "messages_messagebody" VALUES(104,'<a href="/u/3/">Fabio</a> wrote on

<a href="/p/4/">Finding Common Motifs In Sequences</a> : 
I have a few hundred yeast sequences (20-80bp long) and I want to find common motifs (conserved bases at certain indices) in them. I am using a Mac

',3,'Finding Common Motifs In Sequences',NULL,'2009-09-30 22:09:06.677000');
INSERT INTO "messages_messagebody" VALUES(105,'<a href="/u/2/">Istvan Albert</a> wrote on

<a href="/p/5/">Recommend Easy To Use Microarray Clustering Software</a> : 
Feel free to post your favorite clustering tool.

',2,'Recommend Easy To Use Microarray Clustering Software',NULL,'2009-09-30 22:44:22.647000');
INSERT INTO "messages_messagebody" VALUES(106,'<a href="/u/5/">Zhenhai Zhang</a> wrote on

<a href="/p/6/">(Deleted) Test By Zhenhai</a> : 
Hi, I just created my user id a few minutes ago. 

Post this question to see how it works.

',5,'Test By Zhenhai',NULL,'2009-10-01 00:49:39.563000');
INSERT INTO "messages_messagebody" VALUES(107,'<a href="/u/5/">Zhenhai Zhang</a> wrote on

<a href="/p/4/#7">A: Finding Common Motifs In Sequences</a> : 
try this out?

http://fraenkel.mit.edu/webmotifs/form.html

',5,'A: Finding Common Motifs In Sequences',NULL,'2009-10-01 00:55:02.457000');
INSERT INTO "messages_messagebody" VALUES(108,'<a href="/u/6/">Tom Koerber</a> wrote on

<a href="/p/4/#8">A: Finding Common Motifs In Sequences</a> : 
You can also use MEME:  http://meme.sdsc.edu/.

',6,'A: Finding Common Motifs In Sequences',NULL,'2009-10-01 01:32:29.097000');
INSERT INTO "messages_messagebody" VALUES(109,'<a href="/u/7/">Suk211</a> wrote on

<a href="/p/4/#9">A: Finding Common Motifs In Sequences</a> : 

ACGGGCCCGACGATGCGTCGTA

ACGTACGTCGAACCGTCGTCGT

ACGTGCGTCGAAACGTCAGTCG

ACGGGTTCGATCGTCGTCGTCG


may be in Python I will break down the first sequence of required motif length into a sliding window and will search for those list of motifs in the rest of sequences using regular expression in python 
',7,'A: Finding Common Motifs In Sequences',NULL,'2009-10-01 01:35:28.020000');
INSERT INTO "messages_messagebody" VALUES(110,'<a href="/u/10/">Question Bot</a> wrote on

<a href="/p/10/">How To Generate Multi-Nucleotide Occupancy Counts For Each Coordinate Of My Reads?</a> : 
I need to generate nucleotide occupancy counts for each position of a given sequence then summed over each of the input sequences. An example desired output (for di-nucleotide AT):



',10,'How To Generate Multi-Nucleotide Occupancy Counts For Each Coordinate Of My Reads?',NULL,'2009-10-05 21:51:37.043000');
INSERT INTO "messages_messagebody" VALUES(111,'<a href="/u/2/">Istvan Albert</a> wrote on

<a href="/p/10/#11">A: How To Generate Multi-Nucleotide Occupancy Counts For Each Coordinate Of My Read</a> : 
The code snippet below will populate the store dictionary keyed by the nucleotide patterns and values as lists that contain the occupancy for each index. (Updated answer now includes arbitrary lenght nucleotide counts)::

from itertools import count

def pattern_update(sequence, width=2, store={}):

',2,'A: How To Generate Multi-Nucleotide Occupancy Counts For Each Coordinate Of My Read',NULL,'2009-10-05 22:03:01.060000');
INSERT INTO "messages_messagebody" VALUES(112,'<a href="/u/2/">Istvan Albert</a> wrote on

<a href="/p/5/#12">A: Recommend Easy To Use Microarray Clustering Software</a> : 
One of my favorites is the MEV micro-array data analysis tool.
It is simple to use and it has a very large number of features. 

Works well for any type of data. You can also load into it data from a file that is in a simple text format:


GENE1, value1, value2
GENE2, value1, value2


Feel free to p
',2,'A: Recommend Easy To Use Microarray Clustering Software',NULL,'2009-10-05 22:09:57.673000');
INSERT INTO "messages_messagebody" VALUES(113,'<a href="/u/12/">Yang Yang</a> wrote on

<a href="/p/13/">Chip Dna Deep Sequence</a> : 
Hi, everyone,
I am posting this question for my friend.
He is analyzing his CHIP DNA solid deep sequence data, and find out that near 80% reads can not be mapped to the human genome. We are wondering if this high percentage unmapped reads is normal in CHIP DNA deep sequence or there may be something
',12,'Chip Dna Deep Sequence',NULL,'2009-10-07 00:58:10.227000');
INSERT INTO "messages_messagebody" VALUES(114,'<a href="/u/2/">Istvan Albert</a> wrote on

<a href="/p/13/#14">A: Chip Dna Deep Sequence</a> : 
I recall that our first samples that we ran on the Solid sequencer have had bad performance. Not quite an 80% loss but around 40%-60% reads were unmappable (yeast). Some other lab members will hopefully chime in with more details. 

',2,'A: Chip Dna Deep Sequence',NULL,'2009-10-07 02:10:32.730000');
INSERT INTO "messages_messagebody" VALUES(115,'<a href="/u/5/">Zhenhai Zhang</a> wrote on

<a href="/p/13/#15">A: Chip Dna Deep Sequence</a> : 
Hi there,

We have done numbers of SOLiD sequencing run on yeast samples. Normally there are only 30-40 percent of total tags can be uniquely mapped back to yeast genome. 

What I would recommend is do it on solexa. You get much higher quality tags.

cheers,

',5,'A: Chip Dna Deep Sequence',NULL,'2009-10-07 02:41:24.877000');
INSERT INTO "messages_messagebody" VALUES(116,'<a href="/u/13/">Gue Su Chang</a> wrote on

<a href="/p/13/#16">A: Chip Dna Deep Sequence</a> : 
Your 20% mapping yield looks like low for normal ChIP experiment, even for human. Several factors can reduce this mapping yield. I am wondering which kind of ChIP was used in your case. That is, which kind of proteins was ChIPed?

',13,'A: Chip Dna Deep Sequence',NULL,'2009-10-07 04:04:41.747000');
INSERT INTO "messages_messagebody" VALUES(117,'<a href="/u/10/">Question Bot</a> wrote on

<a href="/p/1/#17">A: Site Use Guidelines</a> : 
If you are shy about asking the question on your own behalf submit it to to the Question Bot and it will be posted anonymously. Send email to the Question Bot link at the bottom.

',10,'A: Site Use Guidelines',NULL,'2009-10-07 22:41:44.823000');
INSERT INTO "messages_messagebody" VALUES(118,'<a href="/u/14/">Zhaorong</a> wrote on

<a href="/p/1/#18">A: Site Use Guidelines</a> : 
Hi,

I don&#39;t think a new user can vote on a question or an answer.
The site says I need 15 reputation...

',14,'A: Site Use Guidelines',NULL,'2009-10-09 09:28:20.413000');
INSERT INTO "messages_messagebody" VALUES(119,'<a href="/u/5/">Zhenhai Zhang</a> wrote on

<a href="/p/5/#19">A: Recommend Easy To Use Microarray Clustering Software</a> : 
I would recommend a combination of cluster and treeview.

pretty powerful!

',5,'A: Recommend Easy To Use Microarray Clustering Software',NULL,'2009-10-16 18:25:38.237000');
INSERT INTO "messages_messagebody" VALUES(120,'<a href="/u/15/">Nickey</a> wrote on

<a href="/p/20/">(Deleted) Do You Have To Be A Guy To Dress Up As Boy George</a> : 
any ideas im a girl

',15,'Do You Have To Be A Guy To Dress Up As Boy George',NULL,'2009-10-18 09:22:53.980000');
INSERT INTO "messages_messagebody" VALUES(121,'<a href="/u/15/">Nickey</a> wrote on

<a href="/p/21/">(Deleted) Do You Have To Be A Guy To Dress Up As Boy George</a> : 
any ideas im a girl

',15,'Do You Have To Be A Guy To Dress Up As Boy George',NULL,'2009-10-18 09:23:34.373000');
INSERT INTO "messages_messagebody" VALUES(122,'<a href="/u/16/">Renee</a> wrote on

<a href="/p/22/">Gene Id Conversion Tool</a> : 
Hey,

I was using DAVID (http://david.abcc.ncifcrf.gov/conversion.jsp) to do the gene ID conversion, e.g.conversion between Agilent ID, Genebank accession id and Entrez gene ID, but I found the DAVID database is not updated. Does anyone know a better updataed conversion tool to do this job? Thanks! 
',16,'Gene Id Conversion Tool',NULL,'2009-10-23 23:42:24.427000');
INSERT INTO "messages_messagebody" VALUES(123,'<a href="/u/2/">Istvan Albert</a> wrote on

<a href="/p/22/#23">A: Gene Id Conversion Tool</a> : 
I don&#39;t know of a direct solution myself, but this is a topic that may be of interest for the biological data analysis class that I am teaching. 

If you specify the organism/genomic builds that you are interested in we may be able to generate a full translation list as an in class example or a home
',2,'A: Gene Id Conversion Tool',NULL,'2009-10-24 01:46:45.407000');
INSERT INTO "messages_messagebody" VALUES(124,'<a href="/u/14/">Zhaorong</a> wrote on

<a href="/p/24/">How To Set Shrimp Parameters For Best Sensitivity With 35Bp Colorspace Data?</a> : 
Hi,

I have 35bp Solid colorspace sequencing data, and the actual sequences to be mapped are 20-25bp after removing the linker sequence.

I hope to find all the hits allowing no more than n mismatches (say n=3), not only the best hit.

I know there is a -M option to specify -M sensitivity,35bp. I wo
',14,'How To Set Shrimp Parameters For Best Sensitivity With 35Bp Colorspace Data?',NULL,'2009-12-01 13:13:53.637000');
INSERT INTO "messages_messagebody" VALUES(125,'<a href="/u/18/">Gue Su</a> wrote on

<a href="/p/24/#25">A: How To Set Shrimp Parameters For Best Sensitivity With 35Bp Colorspace Data?</a> : 
I just read the SHRiMP manual again, but I think that their explanation about -M option may not be enough to answer your question. I usually use the &quot;seed&quot; mode by using -s, -n, and -w and the option -M is a new feature of the version 1.3.1, which I have never tried before.

I recommend for you to u
',18,'A: How To Set Shrimp Parameters For Best Sensitivity With 35Bp Colorspace Data?',NULL,'2009-12-01 20:57:35.300000');
INSERT INTO "messages_messagebody" VALUES(126,'<a href="/u/2/">Istvan Albert</a> wrote on

<a href="/p/24/#26">A: How To Set Shrimp Parameters For Best Sensitivity With 35Bp Colorspace Data?</a> : 

  Since my reads are only 20-25bp long,
  should I changed the default 4 spaced
  seeds to 3?


while the shrimp manual says:


We recommend using the default 4 seeds of weight 12 in most cases.


you could try running on a smaller sample and see what happens. 

',2,'A: How To Set Shrimp Parameters For Best Sensitivity With 35Bp Colorspace Data?',NULL,'2009-12-02 00:00:53.443000');
INSERT INTO "messages_messagebody" VALUES(127,'<a href="/u/19/">Mohammed Islaih</a> wrote on

<a href="/p/22/#27">A: Gene Id Conversion Tool</a> : 
The following link has a list of ID conversion tools:

http://hum-molgen.org/NewsGen/08-2009/000020.html

',19,'A: Gene Id Conversion Tool',NULL,'2009-12-09 03:45:54.547000');
INSERT INTO "messages_messagebody" VALUES(128,'<a href="/u/20/">Alex Reynolds</a> wrote on

<a href="/p/28/">Tips On Compiling And Using Meme 4.3 With A Sun Grid Engine Computation Cluster</a> : 
Has anyone compiled and used MEME 4.x for use in a parallel computation environment, based upon operation with a Sun Grid Engine (SGE) cluster?

I can compile the suite and its tests pass. However, when I attempt to use the -p n option, to specify n computation nodes, I get several error messages:


',20,'Tips On Compiling And Using Meme 4.3 With A Sun Grid Engine Computation Cluster',NULL,'2010-01-13 18:59:22.603000');
INSERT INTO "messages_messagebody" VALUES(129,'<a href="/u/2/">Istvan Albert</a> wrote on

<a href="/p/28/#29">A: Tips On Compiling And Using Meme 4.3 With A Sun Grid Engine Computation Cluster</a> : 
This may not be overly useful but it very much sounds like a configuration problem.

Usually there is  configure flag that needs to be set to point to the libraries, something like:

--with-mpidir=MPIDIR
--with-mpicc=MPICC


It also appears that the MEME suite does not support Open MPI (as per insta
',2,'A: Tips On Compiling And Using Meme 4.3 With A Sun Grid Engine Computation Cluster',NULL,'2010-01-14 03:17:50.023000');
INSERT INTO "messages_messagebody" VALUES(130,'<a href="/u/20/">Alex Reynolds</a> wrote on

<a href="/p/2/#30">A: How Do I Convert From Bed Format To Gff Format?</a> : 
Here&#39;s a Perl script I wrote if you wanted to do something local. 

There&#39;s some code in there for translating yeast chromosome names that can be removed, if not needed. I also used a Site feature in the GFF file as the region ID, which might also need tweaking, depending on what features you&#39;re int
',20,'A: How Do I Convert From Bed Format To Gff Format?',NULL,'2010-01-15 14:48:14.660000');
INSERT INTO "messages_messagebody" VALUES(131,'<a href="/u/4/">Jason</a> wrote on

<a href="/p/31/">How Do I Map, Align, And Plot My Solid Results?</a> : 
Hi, I recently performed an RNA immunoprecipitation followed by SOLiD sequencing (50 bp fragmented reads). I haven&#39;t received my first SOLiD sequencing results yet, but I was told I should have them soon. I&#39;ve tried doing my own research on how to map, align, and plot my results but I don&#39;t have a c
',4,'How Do I Map, Align, And Plot My Solid Results?',NULL,'2010-01-22 09:14:17.380000');
INSERT INTO "messages_messagebody" VALUES(132,'<a href="/u/2/">Istvan Albert</a> wrote on

<a href="/p/31/#32">A: How Do I Map, Align, And Plot My Solid Results?</a> : 
Personally I would advise that if you know someone who can partially perform the task you should have them do it, and ask them to explain and show it to you how they&#39;ve done it.

The task at hand is complex. The solution always depends immensely on the particulars of the problem, moreover you will b
',2,'A: How Do I Map, Align, And Plot My Solid Results?',NULL,'2010-01-22 21:13:42.790000');
INSERT INTO "messages_messagebody" VALUES(133,'<a href="/u/23/">Giovanni M Dall&#39;Olio</a> wrote on

<a href="/p/33/">Which Operating System Do You Prefer For Bioinformatics?</a> : 
So, you will probably hate me for asking this question here, as there are lot of forum and blog posts on internet about it and it is also a very subjective question.

However, it may be a starting point for a good discussion, if we don&#39;t flame... Which operating system do you usually use for your wo
',23,'Which Operating System Do You Prefer For Bioinformatics?',NULL,'2010-01-26 21:14:38.623000');
INSERT INTO "messages_messagebody" VALUES(134,'<a href="/u/23/">Giovanni M Dall&#39;Olio</a> wrote on

<a href="/p/34/">Which Are The Best Programming Languages To Study For A Bioinformatician?</a> : 
This is also a very classic question: Which is your favorite programming language in bioinformatics? Which languages would you recommend to a student wishing to enter the world of bioinformatics?

This topic has already been discussed on the Internet, but I think it would be nice to discuss it here.
',23,'Which Are The Best Programming Languages To Study For A Bioinformatician?',NULL,'2010-01-26 21:45:06.737000');
INSERT INTO "messages_messagebody" VALUES(135,'<a href="/u/2/">Istvan Albert</a> wrote on

<a href="/p/33/#35">A: Which Operating System Do You Prefer For Bioinformatics?</a> : 
Often people are limited to their choices by factors outside of their control. One lab that I work with requires the use of Mac computers another is using Windows mostly. Large scale computations seem to be best suited for Linux systems.

Luckily there is a migration towards unified capabilities acr
',2,'A: Which Operating System Do You Prefer For Bioinformatics?',NULL,'2010-01-27 02:23:45.583000');
INSERT INTO "messages_messagebody" VALUES(136,'<a href="/u/2/">Istvan Albert</a> wrote on

<a href="/p/34/#36">A: Which Are The Best Programming Languages To Study For A Bioinformatician?</a> : 
It is important to be considerate and not characterize one particular approach negatively. My favorite quote is:

Programming is pure thought.

Hopefully everyone is able to pick an approach that matches their individual way of thinking. While I myself do not program in Perl, I consider it to be one
',2,'A: Which Are The Best Programming Languages To Study For A Bioinformatician?',NULL,'2010-01-27 02:32:29.450000');
INSERT INTO "messages_messagebody" VALUES(137,'<a href="/u/2/">Istvan Albert</a> wrote on

<a href="/p/33/#37">A: Which Operating System Do You Prefer For Bioinformatics?</a> : 
Tips for installing software om Max OS X:


install the Apple developer tools called Xcode http://developer.apple.com/tools/xcode/
install MacPorts from http://www.macports.org/


You can now easily install everything from command line using the port command. List all available software

port list


',2,'A: Which Operating System Do You Prefer For Bioinformatics?',NULL,'2010-01-27 02:38:24.477000');
INSERT INTO "messages_messagebody" VALUES(138,'<a href="/u/25/">Fabio</a> wrote on

<a href="/p/34/#38">A: Which Are The Best Programming Languages To Study For A Bioinformatician?</a> : 
Any programming language is good as long you know what you&#39;re doing.

',25,'A: Which Are The Best Programming Languages To Study For A Bioinformatician?',NULL,'2010-01-27 05:42:14.167000');
INSERT INTO "messages_messagebody" VALUES(139,'<a href="/u/23/">Giovanni M Dall&#39;Olio</a> wrote on

<a href="/p/22/#39">A: Gene Id Conversion Tool</a> : 
You can also do it with the following services:


uniprot - Click on &#39;Id Mapping&#39; from the home page.
biomart - choose a database and a version, then put the ids you want to convert under Filters-&amp;gt;Id List limit (select the proper input id in the menu), and then the output ids under &#39;Attributes&#39;. 
',23,'A: Gene Id Conversion Tool',NULL,'2010-01-27 16:08:52.870000');
INSERT INTO "messages_messagebody" VALUES(140,'<a href="/u/23/">Giovanni M Dall&#39;Olio</a> wrote on

<a href="/p/4/#40">A: Finding Common Motifs In Sequences</a> : 
Meme has been the first program to be published for doing that.
As an alternative you can find one of the EMBOSS tools; if you are scared by a terminal and want to do it from a web-based interface, you can use the EMBOSS tools from galaxy

',23,'A: Finding Common Motifs In Sequences',NULL,'2010-01-28 21:31:50.470000');
INSERT INTO "messages_messagebody" VALUES(141,'<a href="/u/23/">Giovanni M Dall&#39;Olio</a> wrote on

<a href="/p/41/">How Much Do You Trust Geneontology?</a> : 
GeneOntology is a nice project to provide a standard terminology for genes and gene functions, to help avoid the use of synonyms and wrong spelling when describing a gene.

I have been using the GeneOntology for a while, but honestly I think that it contains many errors and that many terms have not 
',23,'How Much Do You Trust Geneontology?',NULL,'2010-01-28 22:17:37.437000');
INSERT INTO "messages_messagebody" VALUES(142,'<a href="/u/23/">Giovanni M Dall&#39;Olio</a> wrote on

<a href="/p/1/#42">A: Site Use Guidelines</a> : 
The StackExchange websites have been designed for making questions related to programming and technical issues.

For example, for this reason, if you try to write a question which starts with &#39;What is your favorite experience...&#39; you get a disclaimer saying that &#39;your question seems to be probably s
',23,'A: Site Use Guidelines',NULL,'2010-01-28 22:22:04.203000');
INSERT INTO "messages_messagebody" VALUES(143,'<a href="/u/23/">Giovanni M Dall&#39;Olio</a> wrote on

<a href="/p/34/#43">A: Which Are The Best Programming Languages To Study For A Bioinformatician?</a> : 
The choice of a programming language is purely subjective, but when a student asks you which programming language he should start with, you have to make an answer, or at least provide some informations.

I think that a bioinformatician who studies R and at least two or three libraries (lattice/ggplo
',23,'A: Which Are The Best Programming Languages To Study For A Bioinformatician?',NULL,'2010-01-28 23:58:20.500000');
INSERT INTO "messages_messagebody" VALUES(144,'<a href="/u/2/">Istvan Albert</a> wrote on

<a href="/p/41/#44">A: How Much Do You Trust Geneontology?</a> : 
The GO terms and classifications are primarily an based on opinions and a human interpretation of a small group of people of what the current state of the knowledge is.Thus  are more subjective than say experimental measurements would be. 

In fact it is surprising that it works at all; and it does 
',2,'A: How Much Do You Trust Geneontology?',NULL,'2010-01-29 04:41:29.170000');
INSERT INTO "messages_messagebody" VALUES(145,'<a href="/u/4/">Jason</a> wrote on

<a href="/p/41/#45">A: How Much Do You Trust Geneontology?</a> : 
In my experience it&#39;s case by case. In other words just because you are getting significant p-values, does not mean the results are biologically significant. I once submitted clusters of microarray data and received a bunch of hits that were significant by p-value, but really didn&#39;t have a theme. Th
',4,'A: How Much Do You Trust Geneontology?',NULL,'2010-01-29 11:50:26.020000');
INSERT INTO "messages_messagebody" VALUES(146,'<a href="/u/23/">Giovanni M Dall&#39;Olio</a> wrote on

<a href="/p/46/">What Is Your Experience With The String (Interactions) Database?</a> : 
STRING is a database of predicted protein-protein interactions at EMBL. It cluster the results from many sources of protein-protein interactions databases, like Mint, etc.., and it also use the informations from KEGG-pathways and reactome, to provide the best annotations for the interactions of a pr
',23,'What Is Your Experience With The String (Interactions) Database?',NULL,'2010-01-29 17:42:23.043000');
INSERT INTO "messages_messagebody" VALUES(147,'<a href="/u/2/">Istvan Albert</a> wrote on

<a href="/p/46/#47">A: What Is Your Experience With The String (Interactions) Database?</a> : 
I have not used STRING in particular but I have worked with protein interactions before (DIP dataset). I recall that even experimentally produced protein-protein interactions may have very large false positive ratios  (as for false negatives, who knows?) Some papers claim that up to 50% of the inter
',2,'A: What Is Your Experience With The String (Interactions) Database?',NULL,'2010-01-29 20:06:06.180000');
INSERT INTO "messages_messagebody" VALUES(148,'<a href="/u/23/">Giovanni M Dall&#39;Olio</a> wrote on

<a href="/p/48/">Where Can I Get The Secondary Structure Of A Protein?</a> : 
As in the title... I have a protein and I would like to know its secundary structure.
I couldn&#39;t find it in uniprot, althought I tought they had annotations for it there.
In the end I have used a predictor (jpred) but there it should be a database somewhere.

',23,'Where Can I Get The Secondary Structure Of A Protein?',NULL,'2010-02-12 23:33:06.820000');
INSERT INTO "messages_messagebody" VALUES(149,'<a href="/u/2/">Istvan Albert</a> wrote on

<a href="/p/48/#49">A: Where Can I Get The Secondary Structure Of A Protein?</a> : 
Protein structure prediction is a complex issue that is likely to require multiple approaches. There are many methods/tools listed at the 


Expert Protein Analysis System website


',2,'A: Where Can I Get The Secondary Structure Of A Protein?',NULL,'2010-02-13 02:01:49.187000');
INSERT INTO "messages_messagebody" VALUES(150,'<a href="/u/26/">Nicojo</a> wrote on

<a href="/p/48/#50">A: Where Can I Get The Secondary Structure Of A Protein?</a> : 
I think you found the best answer yourself: use a predictor! There are several out there...

You suggest that there should be a Secondary Structure Database. I&#39;m not sure that makes much sense, let me explain my point of view (which may not be that of everyone): most often, the data that is found in
',26,'A: Where Can I Get The Secondary Structure Of A Protein?',NULL,'2010-02-13 03:57:06.680000');
INSERT INTO "messages_messagebody" VALUES(151,'<a href="/u/14/">Zhaorong</a> wrote on

<a href="/p/51/">Turn Off Blast Search On Reverse Complement Strand In Blastn</a> : 
I have a quick question:
How can I turn off search on reverse complement strand of my query nucleotide sequence in blastn?

For example, I don&#39;t want &#39;GUAAAGCCAAAUCUUCGGUUA&#39; to be a hit when I use &#39;UAACCGAAGAUUUGGCUUUAC&#39; as the query.

Maybe I missed it when I read the man page, but I really appreci
',14,'Turn Off Blast Search On Reverse Complement Strand In Blastn',NULL,'2010-02-16 07:13:06.543000');
INSERT INTO "messages_messagebody" VALUES(152,'<a href="/u/2/">Istvan Albert</a> wrote on

<a href="/p/51/#52">A: Turn Off Blast Search On Reverse Complement Strand In Blastn</a> : 
The -S flag can select the strands:

-S  Query strands to search against database 
    (for blast[nx], and tblastx) 3 is both, 1 is top, 2 is bottom [Integer]


',2,'A: Turn Off Blast Search On Reverse Complement Strand In Blastn',NULL,'2010-02-16 08:31:37.060000');
INSERT INTO "messages_messagebody" VALUES(153,'<a href="/u/27/">Allen Yu</a> wrote on

<a href="/p/53/">How To Do Quality Trimming Of Solid Reads In Colour Space?</a> : 
The reads returned from the Solid sequencing provider are littered with dots and some bases have a negative quality value. Does anyone know if there is a good method to extract high quality regions from the reads without distorting the reading of bases in colour space?

',27,'How To Do Quality Trimming Of Solid Reads In Colour Space?',NULL,'2010-02-19 13:29:18.733000');
INSERT INTO "messages_messagebody" VALUES(154,'<a href="/u/27/">Allen Yu</a> wrote on

<a href="/p/31/#54">A: How Do I Map, Align, And Plot My Solid Results?</a> : 
You can try BWA as well:
http://maq.sourceforge.net/bwa-man.shtml

',27,'A: How Do I Map, Align, And Plot My Solid Results?',NULL,'2010-02-19 13:31:24.643000');
INSERT INTO "messages_messagebody" VALUES(155,'<a href="/u/28/">Curious George</a> wrote on

<a href="/p/53/#55">A: How To Do Quality Trimming Of Solid Reads In Colour Space?</a> : 
The Solid Accuracy Enhancer Tool might be useful for this.

',28,'A: How To Do Quality Trimming Of Solid Reads In Colour Space?',NULL,'2010-02-19 19:33:45.533000');
INSERT INTO "messages_messagebody" VALUES(156,'<a href="/u/23/">Giovanni M Dall&#39;Olio</a> wrote on

<a href="/p/56/">How To Get The Sequence Of A Genomic Region From Ucsc?</a> : 
Let&#39;s say I want to download the fasta sequence of the region chr1:100000..200000 from the UCSC browser.
How do you do that? I can&#39;t find a button to &#39;export to fasta&#39; in the UCSC genome browser. I think that the solution is to click on one of the tracks displayed, but I am not sure of which.
If I g
',23,'How To Get The Sequence Of A Genomic Region From Ucsc?',NULL,'2010-02-21 22:13:39.320000');
INSERT INTO "messages_messagebody" VALUES(157,'<a href="/u/2/">Istvan Albert</a> wrote on

<a href="/p/56/#57">A: How To Get The Sequence Of A Genomic Region From Ucsc?</a> : 
The Genome Browser is for visualization.

To get data in many formats use the UCSC Table Browser then select the output format of your choice.

You may also need to select the right group and track to get the data you want.

',2,'A: How To Get The Sequence Of A Genomic Region From Ucsc?',NULL,'2010-02-22 01:11:20.450000');
INSERT INTO "messages_messagebody" VALUES(158,'<a href="/u/23/">Giovanni M Dall&#39;Olio</a> wrote on

<a href="/p/58/">What Is The Best Way To Share Scripts Between Members Of A Lab?</a> : 
One of the most awful problems in my group is avoiding to rewrite scripts that have been already written by others. Since we have different projects and we work with different data, everybody ends up writing its own scripts in his favorite programming language, and it is very frequent to waste an af
',23,'What Is The Best Way To Share Scripts Between Members Of A Lab?',NULL,'2010-02-25 21:39:15.467000');
INSERT INTO "messages_messagebody" VALUES(159,'<a href="/u/30/">Pierre Lindenbaum</a> wrote on

<a href="/p/56/#59">A: How To Get The Sequence Of A Genomic Region From Ucsc?</a> : 
Use the DAS server:

http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr1:100000,200000

',30,'A: How To Get The Sequence Of A Genomic Region From Ucsc?',NULL,'2010-02-25 23:32:53.740000');
INSERT INTO "messages_messagebody" VALUES(160,'<a href="/u/31/">Marcos De Carvalho</a> wrote on

<a href="/p/4/#60">A: Finding Common Motifs In Sequences</a> : 
Some time ago I used SOMBRERO (http://bioinf.nuigalway.ie/sombrero/download.html) with a good degree of success on finding motifs in a very diverse set of sequences. They have a Mac version for download as well as parallel versions for Irix and Linux.

',31,'A: Finding Common Motifs In Sequences',NULL,'2010-02-25 23:51:28.290000');
INSERT INTO "messages_messagebody" VALUES(161,'<a href="/u/31/">Marcos De Carvalho</a> wrote on

<a href="/p/58/#61">A: What Is The Best Way To Share Scripts Between Members Of A Lab?</a> : 
I would recommend you to setup a wiki for your group. If you do not have a server readily you can always use one of the many wiki services available for free like Wikispaces (www.wikispaces.com).

',31,'A: What Is The Best Way To Share Scripts Between Members Of A Lab?',NULL,'2010-02-25 23:54:49.080000');
INSERT INTO "messages_messagebody" VALUES(162,'<a href="/u/2/">Istvan Albert</a> wrote on

<a href="/p/58/#62">A: What Is The Best Way To Share Scripts Between Members Of A Lab?</a> : 
Integrating with the source code management tool is essential, that way when code gets changed everyone can easily get the updated version. Wikis are also a good idea.

',2,'A: What Is The Best Way To Share Scripts Between Members Of A Lab?',NULL,'2010-02-26 00:24:46.380000');
INSERT INTO "messages_messagebody" VALUES(163,'<a href="/u/35/">David Nusinow</a> wrote on

<a href="/p/46/#63">A: What Is Your Experience With The String (Interactions) Database?</a> : 
I&#39;ve been using STRING extensively, but not for protein-protein interactions work. STRING, as you note, is a bit of a mutt in terms of the different data sources it mines. Some that you&#39;re missing include a broad literature-based search, as well as gene expression data sets. So if you&#39;re interested 
',35,'A: What Is Your Experience With The String (Interactions) Database?',NULL,'2010-02-26 03:33:37.903000');
INSERT INTO "messages_messagebody" VALUES(164,'<a href="/u/35/">David Nusinow</a> wrote on

<a href="/p/34/#64">A: Which Are The Best Programming Languages To Study For A Bioinformatician?</a> : 
Perl can be quite lovely if you choose to write it well. If you find yourself in need of writing some perl, I&#39;d highly recommend getting the Perl Best Practices book and going through it to learn how to make your perl code not suck. Essential tools for helping with that are perlcritic and perltidy, 
',35,'A: Which Are The Best Programming Languages To Study For A Bioinformatician?',NULL,'2010-02-26 03:41:22.310000');
INSERT INTO "messages_messagebody" VALUES(165,'<a href="/u/33/">User 6996</a> wrote on

<a href="/p/33/#65">A: Which Operating System Do You Prefer For Bioinformatics?</a> : 
My tip: install Cygwin if you are using Windows 

',33,'A: Which Operating System Do You Prefer For Bioinformatics?',NULL,'2010-02-26 03:50:25.157000');
INSERT INTO "messages_messagebody" VALUES(166,'<a href="/u/33/">User 6996</a> wrote on

<a href="/p/1/#66">A: Site Use Guidelines</a> : 
Who is running this site?

',33,'A: Site Use Guidelines',NULL,'2010-02-26 03:52:10.690000');
INSERT INTO "messages_messagebody" VALUES(167,'<a href="/u/38/">Simon Cockell</a> wrote on

<a href="/p/33/#67">A: Which Operating System Do You Prefer For Bioinformatics?</a> : 
All of the 3 major platforms have their advantages, and I use all 3 practically every day. Mac OS X is my primary desktop OS, for a number of reasons, but mostly because I just seem more productive using it than any of the alternatives. All of my coding work is done over SSH on Linux (almost exclusi
',38,'A: Which Operating System Do You Prefer For Bioinformatics?',NULL,'2010-02-26 15:51:33.997000');
INSERT INTO "messages_messagebody" VALUES(168,'<a href="/u/38/">Simon Cockell</a> wrote on

<a href="/p/58/#68">A: What Is The Best Way To Share Scripts Between Members Of A Lab?</a> : 
If you want to see the code, but also store associated information, such as expected outputs etc, then a wiki probably is the best choice (we prefer DokuWiki here), although this would involve a lot of manual effort to document each script. 

Use of a site such as GitHub would give you version contr
',38,'A: What Is The Best Way To Share Scripts Between Members Of A Lab?',NULL,'2010-02-26 16:07:14.700000');
INSERT INTO "messages_messagebody" VALUES(169,'<a href="/u/30/">Pierre Lindenbaum</a> wrote on

<a href="/p/69/">Using Hdf5 To Store  Bio-Data</a> : 
Hi all,
has anobody ever used the HDF5 API to store some biological data (genotypes...). I know about this kind of reference (BioHDF...)  but I&#39;m looking for some source code I could browse to understand how I can access data faster.

Pierre

PS: hum, I&#39;m a new user. I&#39;m not allowed to add the follo
',30,'Using Hdf5 To Store  Bio-Data',NULL,'2010-02-26 18:50:17.577000');
INSERT INTO "messages_messagebody" VALUES(170,'<a href="/u/40/">Biorelated</a> wrote on

<a href="/p/58/#70">A: What Is The Best Way To Share Scripts Between Members Of A Lab?</a> : 
You might also want to setup a simple snippets database. Navysnip application by Jason Strutz is easy to install and run if you have ruby and rubyonrails installed.

git clone git://github.com/navyrain/navysnip.git
  cd navysnip
  sudo rake gems:install
  rake db:migrate
  ruby script/server
Then vi
',40,'A: What Is The Best Way To Share Scripts Between Members Of A Lab?',NULL,'2010-02-26 19:15:16.893000');
INSERT INTO "messages_messagebody" VALUES(171,'<a href="/u/42/">Fernando Muñiz</a> wrote on

<a href="/p/69/#71">A: Using Hdf5 To Store  Bio-Data</a> : 
Hello Pierre!

I have been talking with the BioHDF guys and from what they tell me, their work will be centered around a number of command-line APIs, written in C, that will address some areas of usage which for now do not seem to overlap. 

I have seen this example on their site:
http://www.hdfgrou
',42,'A: Using Hdf5 To Store  Bio-Data',NULL,'2010-02-26 19:47:41.697000');
INSERT INTO "messages_messagebody" VALUES(172,'<a href="/u/23/">Giovanni M Dall&#39;Olio</a> wrote on

<a href="/p/69/#72">A: Using Hdf5 To Store  Bio-Data</a> : 
Unfortunately I don&#39;t have any example to shows you yet.
I don&#39;t know how to program in C/C++ so I have been looking at two hdf5 wrappers in python, PyTables and H5PY.

PyTables has a database-like approach in which HDF5 is used as a sort of hierarchical database, in which a column can be a table it
',23,'A: Using Hdf5 To Store  Bio-Data',NULL,'2010-02-26 19:54:03.447000');
INSERT INTO "messages_messagebody" VALUES(173,'<a href="/u/2/">Istvan Albert</a> wrote on

<a href="/p/69/#73">A: Using Hdf5 To Store  Bio-Data</a> : 
In the GeneTrack software we have used HDF to store values for each genomic base. Its main advantage over other storage systems was that it was able to return consecutive values with minimal overhead. 

For example it is extremely fast (ms) in retrieving say 100,000 consecutive values starting with 
',2,'A: Using Hdf5 To Store  Bio-Data',NULL,'2010-02-26 19:56:22.753000');
INSERT INTO "messages_messagebody" VALUES(174,'<a href="/u/42/">Fernando Muñiz</a> wrote on

<a href="/p/69/#74">A: Using Hdf5 To Store  Bio-Data</a> : 
What I do have is a netCDF-3 based Java application that I could show you.
NetCDF-3 is basically the same idea as HDF, but quite more limited as it cannot do compound datatypes among other limitations.

But here&#39;s a small test code example to toy with:

package netCDF;

import java.io.File;
import u
',42,'A: Using Hdf5 To Store  Bio-Data',NULL,'2010-02-26 20:26:41.150000');
INSERT INTO "messages_messagebody" VALUES(175,'<a href="/u/26/">Nicojo</a> wrote on

<a href="/p/1/#75">A: Site Use Guidelines</a> : 
This is an excellent initiative: congratulations and thank you for setting it up!

I guess this site will be all the more useful as there are more contributers... So I guess that good questions for the administrator(s) of this site are: 


Do you have a plan for advertising this site/attracting new 
',26,'A: Site Use Guidelines',NULL,'2010-02-26 20:44:01.490000');
INSERT INTO "messages_messagebody" VALUES(176,'<a href="/u/30/">Pierre Lindenbaum</a> wrote on

<a href="/p/76/">Looking For A &#39;Hello World&quot; Plugin For Taverna.</a> : 
Hi all,
I&#39;d like to create a very simple plugin for Taverna 2.0, something very simple like like implementing a &#39;convertDnaToRna&#39;. There is already some source code that can be found on the net e.g. Egon Willighagen&#39;s code at http://github.com/egonw/cdk-taverna but it requires to know Maven and.... 
',30,'Looking For A ''Hello World" Plugin For Taverna.',NULL,'2010-02-26 21:37:25.810000');
INSERT INTO "messages_messagebody" VALUES(177,'<a href="/u/44/">Andrew</a> wrote on

<a href="/p/77/">How Do I Convert An Illumina Export File To Bed?</a> : 
I have some illumina data generated from the latest version of the illumina pipeline (1.6.0) I need to convert my data into BED to view in ucsc genome browser.

This seems like it should be a fairly common task, however, I am unable to find any scripts to convert my data.

',44,'How Do I Convert An Illumina Export File To Bed?',NULL,'2010-02-26 21:49:18.627000');
INSERT INTO "messages_messagebody" VALUES(178,'<a href="/u/10/">Question Bot</a> wrote on

<a href="/p/77/#78">A: How Do I Convert An Illumina Export File To Bed?</a> : 
I found a script on another site, Uses perl but I have not checked for correctness:

#!/usr/bin/perl

use strict;
use warnings;
use IO::File;

my $filename = shift @ARGV;
die &quot;Usage\n\tperl sorted2bed.pl s_X_sorted.txt &amp;gt; s_X_sorted.bed\n&quot; unless $filename;
chomp $filename;

my $fh = new IO::File;
',10,'A: How Do I Convert An Illumina Export File To Bed?',NULL,'2010-02-26 22:15:18.997000');
INSERT INTO "messages_messagebody" VALUES(179,'<a href="/u/23/">Giovanni M Dall&#39;Olio</a> wrote on

<a href="/p/79/">How To Organize A Pipeline Of Small Scripts Together?</a> : 
In bioinformatics it is very common to end up with a lot of small scripts, each one with a different scope - plotting a chart, converting a file into another format, execute small operations - so it is very important to have a good way to clue them together, to define which should be executed before
',23,'How To Organize A Pipeline Of Small Scripts Together?',NULL,'2010-02-26 22:49:35.150000');
INSERT INTO "messages_messagebody" VALUES(180,'<a href="/u/2/">Istvan Albert</a> wrote on

<a href="/p/79/#80">A: How To Organize A Pipeline Of Small Scripts Together?</a> : 
I don&#39;t have personal experience with this package but it is something that I plan to explore in the near future:

Ruffus  a lightweight python module to run computational pipelines. 

',2,'A: How To Organize A Pipeline Of Small Scripts Together?',NULL,'2010-02-26 22:58:09.550000');
INSERT INTO "messages_messagebody" VALUES(181,'<a href="/u/23/">Giovanni M Dall&#39;Olio</a> wrote on

<a href="/p/79/#81">A: How To Organize A Pipeline Of Small Scripts Together?</a> : 
My favorite way of defining pipelines is by writing Makefiles, about which you can find a very good introduction in Software Carpentry for Bioinformatics: http://swc.scipy.org/lec/build.html .

Although they have been originally developed for compiling programs, Makefiles allow to define which opera
',23,'A: How To Organize A Pipeline Of Small Scripts Together?',NULL,'2010-02-26 23:03:52.940000');
INSERT INTO "messages_messagebody" VALUES(182,'<a href="/u/47/">Chris</a> wrote on

<a href="/p/79/#82">A: How To Organize A Pipeline Of Small Scripts Together?</a> : 
Since I work a lot with Python, I usually write a wrapper method that embeds the external script/program, i.e. calls it, parses its output and returns the desired information. The &#39;glueing&#39; of several such methods then takes place within my Python code that calls all these wrappers. I guess that&#39;s a
',47,'A: How To Organize A Pipeline Of Small Scripts Together?',NULL,'2010-02-26 23:04:53.577000');
INSERT INTO "messages_messagebody" VALUES(183,'<a href="/u/7/">Suk211</a> wrote on

<a href="/p/48/#83">A: Where Can I Get The Secondary Structure Of A Protein?</a> : 
If you have the PDB file then you can use the standard tool called DSSP , it is supposed to be the gold standard for obtaining secondary structure. In case you just have sequence then I personally prefer PSIPRED , it takes evolutionary information into account to predict the secondary structure . Ac
',7,'A: Where Can I Get The Secondary Structure Of A Protein?',NULL,'2010-02-26 23:16:41.967000');
INSERT INTO "messages_messagebody" VALUES(184,'<a href="/u/53/">Michael Barton</a> wrote on

<a href="/p/79/#84">A: How To Organize A Pipeline Of Small Scripts Together?</a> : 
My answer would be: don&#39;t bother. I&#39;ve often found that much of the scripts I write are never used again after the initial use. Therefore spending time using a complex framework that considers dependency between scripts is a waste because the results might be negative and you never visit the analysi
',53,'A: How To Organize A Pipeline Of Small Scripts Together?',NULL,'2010-02-27 02:16:31.970000');
INSERT INTO "messages_messagebody" VALUES(185,'<a href="/u/24/">Etal</a> wrote on

<a href="/p/58/#85">A: What Is The Best Way To Share Scripts Between Members Of A Lab?</a> : 
My lab uses a network-attached storage unit which every Linux workstation mounts by NFS at startup. It was reasonably cheap -- a couple hundred dollars per TB. We also keep copies of public databases on there. We put data sets on there as we&#39;re working on them, and also put the more important script
',24,'A: What Is The Best Way To Share Scripts Between Members Of A Lab?',NULL,'2010-02-27 02:27:32.303000');
INSERT INTO "messages_messagebody" VALUES(186,'<a href="/u/7/">Suk211</a> wrote on

<a href="/p/58/#86">A: What Is The Best Way To Share Scripts Between Members Of A Lab?</a> : 
This might be useful .

A Quick Guide to Organizing Computational Biology Projects

',7,'A: What Is The Best Way To Share Scripts Between Members Of A Lab?',NULL,'2010-02-27 02:44:04.547000');
INSERT INTO "messages_messagebody" VALUES(187,'<a href="/u/55/">Michael Dondrup</a> wrote on

<a href="/p/69/#87">A: Using Hdf5 To Store  Bio-Data</a> : 
There is also a Perl binding to HDF5: PDL::IO::HDF5

http://search.cpan.org/~cerney/PDL-IO-HDF5-0.5/
This requires the Perl Data Language (PDL) package. The way, data-structures can be handled, sub-ranges of data can be defined  an data can be manipulated is actually very elegant in PDL such that co
',55,'A: Using Hdf5 To Store  Bio-Data',NULL,'2010-02-28 03:00:56.140000');
INSERT INTO "messages_messagebody" VALUES(188,'<a href="/u/23/">Giovanni M Dall&#39;Olio</a> wrote on

<a href="/p/88/">Agile Programming For Bioinformaticians - Any Suggestions?</a> : 
I am planning to prepare a talk for my workmates, to introduce them the basics of some agile programming methodology, which I think could give us good ideas to improve our working as a team.

My idea was to take inspiration from extreme programming and explain the rules I like the most: use of A7 ca
',23,'Agile Programming For Bioinformaticians - Any Suggestions?',NULL,'2010-03-01 19:51:12.300000');
INSERT INTO "messages_messagebody" VALUES(189,'<a href="/u/2/">Istvan Albert</a> wrote on

<a href="/p/88/#89">A: Agile Programming For Bioinformaticians - Any Suggestions?</a> : 
I think the approach is unsuited for individuals who are not comfortable with programming in general. There is a long way to go until someone becomes confident in their abilities. Before that this approach is not only ineffective, it might be even be detrimental.

Instead what helps most is transpar
',2,'A: Agile Programming For Bioinformaticians - Any Suggestions?',NULL,'2010-03-01 20:06:55.973000');
INSERT INTO "messages_messagebody" VALUES(190,'<a href="/u/10/">Question Bot</a> wrote on

<a href="/p/90/">Computing The Reverse And Complement Of A Sequence With Biopython</a> : 
An example that computes the reverse complement of a sequence with BioPython

#
# Reverse complement example with BioPython
#

from Bio.Seq import Seq

# a separate function to reverse strings (or other iterables)
def rev(it):
    &quot;Reverses an interable and returns it as a string&quot;
    return &#39;&#39;.join
',10,'Computing The Reverse And Complement Of A Sequence With Biopython',NULL,'2010-03-01 20:26:09.937000');
INSERT INTO "messages_messagebody" VALUES(191,'<a href="/u/39/">Konrad</a> wrote on

<a href="/p/88/#91">A: Agile Programming For Bioinformaticians - Any Suggestions?</a> : 
I would suggest to have a look at Scrum, too. Certain parts would help not only bioinformations. For example estimating the time expenditure of tasks and the resulting burn down charts can be really helpful to see if something is stuck especially when working together on bigger projects.The daily sc
',39,'A: Agile Programming For Bioinformaticians - Any Suggestions?',NULL,'2010-03-01 20:41:28.687000');
INSERT INTO "messages_messagebody" VALUES(192,'<a href="/u/10/">Question Bot</a> wrote on

<a href="/p/92/">(Closed) Computing The Reverse And Complement Of A Sequence With Pygr</a> : 
Computing the reverse complement with the Pygr bioinformatics framework:

#
# Reverse complement example with pygr
#

from pygr.sequence import Sequence

# needs a separate function to reverse strings
def rev(it):
    &quot;Reverses an interable and returns it as a string&quot;
    return &#39;&#39;.join(reversed(it)
',10,'Computing The Reverse And Complement Of A Sequence With Pygr',NULL,'2010-03-01 20:48:25.830000');
INSERT INTO "messages_messagebody" VALUES(193,'<a href="/u/45/">Alex</a> wrote on

<a href="/p/5/#93">A: Recommend Easy To Use Microarray Clustering Software</a> : 
Possibly related:
http://mmc.gnets.ncsu.edu/

',45,'A: Recommend Easy To Use Microarray Clustering Software',NULL,'2010-03-01 20:51:52.190000');
INSERT INTO "messages_messagebody" VALUES(194,'<a href="/u/23/">Giovanni M Dall&#39;Olio</a> wrote on

<a href="/p/90/#94">A: Computing The Reverse And Complement Of A Sequence With Biopython</a> : 
Wouldn&#39;t it better to have a single question titled &#39;How to compute the reverse complement with python&#39; and put all the examples as different answers? Otherwise it seems a bit confusing..

',23,'A: Computing The Reverse And Complement Of A Sequence With Biopython',NULL,'2010-03-01 22:07:15.283000');
INSERT INTO "messages_messagebody" VALUES(195,'<a href="/u/24/">Etal</a> wrote on

<a href="/p/79/#95">A: How To Organize A Pipeline Of Small Scripts Together?</a> : 
The most important thing for me has been keeping a README file at the top of each project directory, where I write down not just how to run the scripts, but why I wrote them in the first place -- coming back to a project after a several-month lull, it&#39;s remarkable difficult to figure out what all th
',24,'A: How To Organize A Pipeline Of Small Scripts Together?',NULL,'2010-03-01 22:32:26.107000');
INSERT INTO "messages_messagebody" VALUES(196,'<a href="/u/24/">Etal</a> wrote on

<a href="/p/90/#96">A: Computing The Reverse And Complement Of A Sequence With Biopython</a> : 
The Bio.Seq module provides two easy ways to get the complement and reverse complement from a sequence:


If you have a string, use the functions complement(dna) and reverse_complement(dna)
If you have a Seq object, use its methods with the same names: dna.complement() and dna.reverse_complement


T
',24,'A: Computing The Reverse And Complement Of A Sequence With Biopython',NULL,'2010-03-01 22:52:22.510000');
INSERT INTO "messages_messagebody" VALUES(197,'<a href="/u/59/">Daniel Swan</a> wrote on

<a href="/p/22/#97">A: Gene Id Conversion Tool</a> : 
http://idconverter.bioinfo.cnio.es/

Is another possible solution to this, although you might find this is not as up to date as you might like either.

',59,'A: Gene Id Conversion Tool',NULL,'2010-03-01 23:05:13.920000');
INSERT INTO "messages_messagebody" VALUES(198,'<a href="/u/55/">Michael Dondrup</a> wrote on

<a href="/p/22/#98">A: Gene Id Conversion Tool</a> : 
BioMart has already been mentioned. It can do much more than ID conversion but it is very useful for conversion purposes, it is regularly updated and you can select different genome builds and all kinds of genomic features. It seems to me that you wish to retrieve GeneIDs linked to Affymetrix IDs. T
',55,'A: Gene Id Conversion Tool',NULL,'2010-03-02 01:51:25.980000');
INSERT INTO "messages_messagebody" VALUES(199,'<a href="/u/10/">Question Bot</a> wrote on

<a href="/p/99/">Fast Interval Intersection Methodologies</a> : 
Most genomic annotations are specified as intervals along the genome. 


Interval trees have been known to provide an efficient datastructure that allows for very fast overlap querying. 
Nested Containment Lists have been proposed as an even faster alternative 


Provide code examples in your progra
',10,'Fast Interval Intersection Methodologies',NULL,'2010-03-02 02:01:20.047000');
INSERT INTO "messages_messagebody" VALUES(200,'<a href="/u/2/">Istvan Albert</a> wrote on

<a href="/p/99/#100">A: Fast Interval Intersection Methodologies</a> : 
This code example generates 10,000 intervals then queries them for overlapping regions. Requires only the presence of Python.

The code below requires the either the installation of the bx python package or alternatively you may just download the quicksect.py module and place it next to the script i
',2,'A: Fast Interval Intersection Methodologies',NULL,'2010-03-02 02:06:18.193000');
INSERT INTO "messages_messagebody" VALUES(201,'<a href="/u/1/">Biostar Community</a> wrote on
<a href="/p/101/">How to find motifs with Galaxy?</a> :
I have a few hundred yeast sequences (20-80bp long) and I want to find common motifs (conserved bases at certain indices) in them. I am using a Mac',1,'How to find motifs with Galaxy?',NULL,'2014-05-09 14:39:20.047000');
INSERT INTO "messages_messagebody" VALUES(202,'Hello <b>john</b>! Welcome to Biostar <i class="fa fa-globe"></i>!

<p style="padding-top:10px;">
    We hope you''ll find the site useful and our community friendly.
</p>
<p>
    We encourage you to participate:
    <b>upvote</b> posts that you like, ask <b>questions</b>, write <b>comments</b> and <b>answers</b>, follow posts,
    bookmark content. You can keep up with the lates content via email (see your user account settings) or RSS feeds.
</p>
<p>
    Together, by sharing the know-how, personal experiences, tips and tricks of the trade
    we can become more productive and successful as scientists.
</p>

',1,'Welcome!',NULL,'2014-05-09 14:39:23.890000');
CREATE TABLE "messages_message" (
    "id" integer NOT NULL PRIMARY KEY,
    "user_id" integer NOT NULL,
    "body_id" integer NOT NULL REFERENCES "messages_messagebody" ("id"),
    "type" integer NOT NULL,
    "unread" bool NOT NULL,
    "sent_at" datetime
);
INSERT INTO "messages_message" VALUES(1,1,1,0,1,'2014-04-29 15:02:11.531000');
INSERT INTO "messages_message" VALUES(2,2,2,0,1,'2014-04-29 15:02:12.464000');
INSERT INTO "messages_message" VALUES(3,3,3,0,1,'2014-04-29 15:02:12.506000');
INSERT INTO "messages_message" VALUES(4,4,4,0,1,'2014-04-29 15:02:12.516000');
INSERT INTO "messages_message" VALUES(5,5,5,0,1,'2014-04-29 15:02:12.528000');
INSERT INTO "messages_message" VALUES(6,6,6,0,1,'2014-04-29 15:02:12.539000');
INSERT INTO "messages_message" VALUES(7,7,7,0,1,'2014-04-29 15:02:12.549000');
INSERT INTO "messages_message" VALUES(8,8,8,0,1,'2014-04-29 15:02:12.562000');
INSERT INTO "messages_message" VALUES(9,9,9,0,1,'2014-04-29 15:02:12.572000');
INSERT INTO "messages_message" VALUES(10,10,10,0,1,'2014-04-29 15:02:12.598000');
INSERT INTO "messages_message" VALUES(11,11,11,0,1,'2014-04-29 15:02:12.608000');
INSERT INTO "messages_message" VALUES(12,12,12,0,1,'2014-04-29 15:02:12.619000');
INSERT INTO "messages_message" VALUES(13,13,13,0,1,'2014-04-29 15:02:12.630000');
INSERT INTO "messages_message" VALUES(14,14,14,0,1,'2014-04-29 15:02:12.640000');
INSERT INTO "messages_message" VALUES(15,15,15,0,1,'2014-04-29 15:02:12.652000');
INSERT INTO "messages_message" VALUES(16,16,16,0,1,'2014-04-29 15:02:12.661000');
INSERT INTO "messages_message" VALUES(17,17,17,0,1,'2014-04-29 15:02:12.672000');
INSERT INTO "messages_message" VALUES(18,18,18,0,1,'2014-04-29 15:02:12.683000');
INSERT INTO "messages_message" VALUES(19,19,19,0,1,'2014-04-29 15:02:12.693000');
INSERT INTO "messages_message" VALUES(20,20,20,0,1,'2014-04-29 15:02:12.704000');
INSERT INTO "messages_message" VALUES(21,21,21,0,1,'2014-04-29 15:02:12.714000');
INSERT INTO "messages_message" VALUES(22,22,22,0,1,'2014-04-29 15:02:12.724000');
INSERT INTO "messages_message" VALUES(23,23,23,0,1,'2014-04-29 15:02:12.735000');
INSERT INTO "messages_message" VALUES(24,24,24,0,1,'2014-04-29 15:02:12.746000');
INSERT INTO "messages_message" VALUES(25,25,25,0,1,'2014-04-29 15:02:12.756000');
INSERT INTO "messages_message" VALUES(26,26,26,0,1,'2014-04-29 15:02:12.766000');
INSERT INTO "messages_message" VALUES(27,27,27,0,1,'2014-04-29 15:02:12.777000');
INSERT INTO "messages_message" VALUES(28,28,28,0,1,'2014-04-29 15:02:12.788000');
INSERT INTO "messages_message" VALUES(29,29,29,0,1,'2014-04-29 15:02:12.799000');
INSERT INTO "messages_message" VALUES(30,30,30,0,1,'2014-04-29 15:02:12.810000');
INSERT INTO "messages_message" VALUES(31,31,31,0,1,'2014-04-29 15:02:12.823000');
INSERT INTO "messages_message" VALUES(32,32,32,0,1,'2014-04-29 15:02:12.836000');
INSERT INTO "messages_message" VALUES(33,33,33,0,1,'2014-04-29 15:02:12.848000');
INSERT INTO "messages_message" VALUES(34,34,34,0,1,'2014-04-29 15:02:12.859000');
INSERT INTO "messages_message" VALUES(35,35,35,0,1,'2014-04-29 15:02:12.869000');
INSERT INTO "messages_message" VALUES(36,36,36,0,1,'2014-04-29 15:02:12.880000');
INSERT INTO "messages_message" VALUES(37,37,37,0,1,'2014-04-29 15:02:12.894000');
INSERT INTO "messages_message" VALUES(38,38,38,0,1,'2014-04-29 15:02:12.905000');
INSERT INTO "messages_message" VALUES(39,39,39,0,1,'2014-04-29 15:02:12.918000');
INSERT INTO "messages_message" VALUES(40,40,40,0,1,'2014-04-29 15:02:12.928000');
INSERT INTO "messages_message" VALUES(41,41,41,0,1,'2014-04-29 15:02:12.939000');
INSERT INTO "messages_message" VALUES(42,42,42,0,1,'2014-04-29 15:02:12.953000');
INSERT INTO "messages_message" VALUES(43,43,43,0,1,'2014-04-29 15:02:12.963000');
INSERT INTO "messages_message" VALUES(44,44,44,0,1,'2014-04-29 15:02:12.974000');
INSERT INTO "messages_message" VALUES(45,45,45,0,1,'2014-04-29 15:02:12.985000');
INSERT INTO "messages_message" VALUES(46,46,46,0,1,'2014-04-29 15:02:12.995000');
INSERT INTO "messages_message" VALUES(47,47,47,0,1,'2014-04-29 15:02:13.006000');
INSERT INTO "messages_message" VALUES(48,48,48,0,1,'2014-04-29 15:02:13.017000');
INSERT INTO "messages_message" VALUES(49,49,49,0,1,'2014-04-29 15:02:13.028000');
INSERT INTO "messages_message" VALUES(50,50,50,0,1,'2014-04-29 15:02:13.039000');
INSERT INTO "messages_message" VALUES(51,51,51,0,1,'2014-04-29 15:02:13.055000');
INSERT INTO "messages_message" VALUES(52,52,52,0,1,'2014-04-29 15:02:13.065000');
INSERT INTO "messages_message" VALUES(53,53,53,0,1,'2014-04-29 15:02:13.077000');
INSERT INTO "messages_message" VALUES(54,54,54,0,1,'2014-04-29 15:02:13.090000');
INSERT INTO "messages_message" VALUES(55,55,55,0,1,'2014-04-29 15:02:13.101000');
INSERT INTO "messages_message" VALUES(56,56,56,0,1,'2014-04-29 15:02:13.114000');
INSERT INTO "messages_message" VALUES(57,57,57,0,1,'2014-04-29 15:02:13.127000');
INSERT INTO "messages_message" VALUES(58,58,58,0,1,'2014-04-29 15:02:13.139000');
INSERT INTO "messages_message" VALUES(59,59,59,0,1,'2014-04-29 15:02:13.150000');
INSERT INTO "messages_message" VALUES(60,60,60,0,1,'2014-04-29 15:02:13.165000');
INSERT INTO "messages_message" VALUES(61,61,61,0,1,'2014-04-29 15:02:13.176000');
INSERT INTO "messages_message" VALUES(62,62,62,0,1,'2014-04-29 15:02:13.189000');
INSERT INTO "messages_message" VALUES(63,63,63,0,1,'2014-04-29 15:02:13.201000');
INSERT INTO "messages_message" VALUES(64,64,64,0,1,'2014-04-29 15:02:13.212000');
INSERT INTO "messages_message" VALUES(65,65,65,0,1,'2014-04-29 15:02:13.239000');
INSERT INTO "messages_message" VALUES(66,66,66,0,1,'2014-04-29 15:02:13.259000');
INSERT INTO "messages_message" VALUES(67,67,67,0,1,'2014-04-29 15:02:13.273000');
INSERT INTO "messages_message" VALUES(68,68,68,0,1,'2014-04-29 15:02:13.288000');
INSERT INTO "messages_message" VALUES(69,69,69,0,1,'2014-04-29 15:02:13.301000');
INSERT INTO "messages_message" VALUES(70,70,70,0,1,'2014-04-29 15:02:13.315000');
INSERT INTO "messages_message" VALUES(71,71,71,0,1,'2014-04-29 15:02:13.339000');
INSERT INTO "messages_message" VALUES(72,72,72,0,1,'2014-04-29 15:02:13.354000');
INSERT INTO "messages_message" VALUES(73,73,73,0,1,'2014-04-29 15:02:13.370000');
INSERT INTO "messages_message" VALUES(74,74,74,0,1,'2014-04-29 15:02:13.385000');
INSERT INTO "messages_message" VALUES(75,75,75,0,1,'2014-04-29 15:02:13.397000');
INSERT INTO "messages_message" VALUES(76,76,76,0,1,'2014-04-29 15:02:13.409000');
INSERT INTO "messages_message" VALUES(77,77,77,0,1,'2014-04-29 15:02:13.420000');
INSERT INTO "messages_message" VALUES(78,78,78,0,1,'2014-04-29 15:02:13.440000');
INSERT INTO "messages_message" VALUES(79,79,79,0,1,'2014-04-29 15:02:13.451000');
INSERT INTO "messages_message" VALUES(80,80,80,0,1,'2014-04-29 15:02:13.463000');
INSERT INTO "messages_message" VALUES(81,81,81,0,1,'2014-04-29 15:02:13.473000');
INSERT INTO "messages_message" VALUES(82,82,82,0,1,'2014-04-29 15:02:13.489000');
INSERT INTO "messages_message" VALUES(83,83,83,0,1,'2014-04-29 15:02:13.500000');
INSERT INTO "messages_message" VALUES(84,84,84,0,1,'2014-04-29 15:02:13.511000');
INSERT INTO "messages_message" VALUES(85,85,85,0,1,'2014-04-29 15:02:13.522000');
INSERT INTO "messages_message" VALUES(86,86,86,0,1,'2014-04-29 15:02:13.535000');
INSERT INTO "messages_message" VALUES(87,87,87,0,1,'2014-04-29 15:02:13.546000');
INSERT INTO "messages_message" VALUES(88,88,88,0,1,'2014-04-29 15:02:13.565000');
INSERT INTO "messages_message" VALUES(89,89,89,0,1,'2014-04-29 15:02:13.580000');
INSERT INTO "messages_message" VALUES(90,90,90,0,1,'2014-04-29 15:02:13.593000');
INSERT INTO "messages_message" VALUES(91,91,91,0,1,'2014-04-29 15:02:13.619000');
INSERT INTO "messages_message" VALUES(92,92,92,0,1,'2014-04-29 15:02:13.637000');
INSERT INTO "messages_message" VALUES(93,93,93,0,1,'2014-04-29 15:02:13.655000');
INSERT INTO "messages_message" VALUES(94,94,94,0,1,'2014-04-29 15:02:13.674000');
INSERT INTO "messages_message" VALUES(95,95,95,0,1,'2014-04-29 15:02:13.688000');
INSERT INTO "messages_message" VALUES(96,96,96,0,1,'2014-04-29 15:02:13.700000');
INSERT INTO "messages_message" VALUES(97,97,97,0,1,'2014-04-29 15:02:13.725000');
INSERT INTO "messages_message" VALUES(98,98,98,0,1,'2014-04-29 15:02:13.737000');
INSERT INTO "messages_message" VALUES(99,99,99,0,1,'2014-04-29 15:02:13.761000');
INSERT INTO "messages_message" VALUES(100,100,100,0,1,'2014-04-29 15:02:13.774000');
INSERT INTO "messages_message" VALUES(101,3,107,0,1,'2009-10-01 00:55:02.457000');
INSERT INTO "messages_message" VALUES(102,3,108,0,1,'2009-10-01 01:32:29.097000');
INSERT INTO "messages_message" VALUES(103,5,108,0,1,'2009-10-01 01:32:29.097000');
INSERT INTO "messages_message" VALUES(104,3,109,0,1,'2009-10-01 01:35:28.020000');
INSERT INTO "messages_message" VALUES(105,5,109,0,1,'2009-10-01 01:35:28.020000');
INSERT INTO "messages_message" VALUES(106,6,109,0,1,'2009-10-01 01:35:28.020000');
INSERT INTO "messages_message" VALUES(107,10,111,0,1,'2009-10-05 22:03:01.060000');
INSERT INTO "messages_message" VALUES(108,12,114,0,1,'2009-10-07 02:10:32.730000');
INSERT INTO "messages_message" VALUES(109,12,115,0,1,'2009-10-07 02:41:24.877000');
INSERT INTO "messages_message" VALUES(110,2,115,0,1,'2009-10-07 02:41:24.877000');
INSERT INTO "messages_message" VALUES(111,12,116,0,1,'2009-10-07 04:04:41.747000');
INSERT INTO "messages_message" VALUES(112,2,116,0,1,'2009-10-07 04:04:41.747000');
INSERT INTO "messages_message" VALUES(113,5,116,0,1,'2009-10-07 04:04:41.747000');
INSERT INTO "messages_message" VALUES(114,2,117,0,1,'2009-10-07 22:41:44.823000');
INSERT INTO "messages_message" VALUES(115,2,118,0,1,'2009-10-09 09:28:20.413000');
INSERT INTO "messages_message" VALUES(116,10,118,0,1,'2009-10-09 09:28:20.413000');
INSERT INTO "messages_message" VALUES(117,2,119,0,1,'2009-10-16 18:25:38.237000');
INSERT INTO "messages_message" VALUES(118,16,123,0,1,'2009-10-24 01:46:45.407000');
INSERT INTO "messages_message" VALUES(119,14,125,0,1,'2009-12-01 20:57:35.300000');
INSERT INTO "messages_message" VALUES(120,14,126,0,1,'2009-12-02 00:00:53.443000');
INSERT INTO "messages_message" VALUES(121,18,126,0,1,'2009-12-02 00:00:53.443000');
INSERT INTO "messages_message" VALUES(122,16,127,0,1,'2009-12-09 03:45:54.547000');
INSERT INTO "messages_message" VALUES(123,2,127,0,1,'2009-12-09 03:45:54.547000');
INSERT INTO "messages_message" VALUES(124,20,129,0,1,'2010-01-14 03:17:50.023000');
INSERT INTO "messages_message" VALUES(125,2,130,0,1,'2010-01-15 14:48:14.660000');
INSERT INTO "messages_message" VALUES(126,4,132,0,1,'2010-01-22 21:13:42.790000');
INSERT INTO "messages_message" VALUES(127,23,135,0,1,'2010-01-27 02:23:45.583000');
INSERT INTO "messages_message" VALUES(128,23,136,0,1,'2010-01-27 02:32:29.450000');
INSERT INTO "messages_message" VALUES(129,23,137,0,1,'2010-01-27 02:38:24.477000');
INSERT INTO "messages_message" VALUES(130,23,138,0,1,'2010-01-27 05:42:14.167000');
INSERT INTO "messages_message" VALUES(131,2,138,0,1,'2010-01-27 05:42:14.167000');
INSERT INTO "messages_message" VALUES(132,16,139,0,1,'2010-01-27 16:08:52.870000');
INSERT INTO "messages_message" VALUES(133,2,139,0,1,'2010-01-27 16:08:52.870000');
INSERT INTO "messages_message" VALUES(134,19,139,0,1,'2010-01-27 16:08:52.870000');
INSERT INTO "messages_message" VALUES(135,3,140,0,1,'2010-01-28 21:31:50.470000');
INSERT INTO "messages_message" VALUES(136,5,140,0,1,'2010-01-28 21:31:50.470000');
INSERT INTO "messages_message" VALUES(137,6,140,0,1,'2010-01-28 21:31:50.470000');
INSERT INTO "messages_message" VALUES(138,7,140,0,1,'2010-01-28 21:31:50.470000');
INSERT INTO "messages_message" VALUES(139,2,142,0,1,'2010-01-28 22:22:04.203000');
INSERT INTO "messages_message" VALUES(140,10,142,0,1,'2010-01-28 22:22:04.203000');
INSERT INTO "messages_message" VALUES(141,14,142,0,1,'2010-01-28 22:22:04.203000');
INSERT INTO "messages_message" VALUES(142,2,143,0,1,'2010-01-28 23:58:20.500000');
INSERT INTO "messages_message" VALUES(143,25,143,0,1,'2010-01-28 23:58:20.500000');
INSERT INTO "messages_message" VALUES(144,23,144,0,1,'2010-01-29 04:41:29.170000');
INSERT INTO "messages_message" VALUES(145,23,145,0,1,'2010-01-29 11:50:26.020000');
INSERT INTO "messages_message" VALUES(146,2,145,0,1,'2010-01-29 11:50:26.020000');
INSERT INTO "messages_message" VALUES(147,23,147,0,1,'2010-01-29 20:06:06.180000');
INSERT INTO "messages_message" VALUES(148,23,149,0,1,'2010-02-13 02:01:49.187000');
INSERT INTO "messages_message" VALUES(149,23,150,0,1,'2010-02-13 03:57:06.680000');
INSERT INTO "messages_message" VALUES(150,2,150,0,1,'2010-02-13 03:57:06.680000');
INSERT INTO "messages_message" VALUES(151,14,152,0,1,'2010-02-16 08:31:37.060000');
INSERT INTO "messages_message" VALUES(152,4,154,0,1,'2010-02-19 13:31:24.643000');
INSERT INTO "messages_message" VALUES(153,2,154,0,1,'2010-02-19 13:31:24.643000');
INSERT INTO "messages_message" VALUES(154,27,155,0,1,'2010-02-19 19:33:45.533000');
INSERT INTO "messages_message" VALUES(155,23,157,0,1,'2010-02-22 01:11:20.450000');
INSERT INTO "messages_message" VALUES(156,23,159,0,1,'2010-02-25 23:32:53.740000');
INSERT INTO "messages_message" VALUES(157,2,159,0,1,'2010-02-25 23:32:53.740000');
INSERT INTO "messages_message" VALUES(158,3,160,0,1,'2010-02-25 23:51:28.290000');
INSERT INTO "messages_message" VALUES(159,5,160,0,1,'2010-02-25 23:51:28.290000');
INSERT INTO "messages_message" VALUES(160,6,160,0,1,'2010-02-25 23:51:28.290000');
INSERT INTO "messages_message" VALUES(161,7,160,0,1,'2010-02-25 23:51:28.290000');
INSERT INTO "messages_message" VALUES(162,23,160,0,1,'2010-02-25 23:51:28.290000');
INSERT INTO "messages_message" VALUES(163,23,161,0,1,'2010-02-25 23:54:49.080000');
INSERT INTO "messages_message" VALUES(164,23,162,0,1,'2010-02-26 00:24:46.380000');
INSERT INTO "messages_message" VALUES(165,31,162,0,1,'2010-02-26 00:24:46.380000');
INSERT INTO "messages_message" VALUES(166,23,163,0,1,'2010-02-26 03:33:37.903000');
INSERT INTO "messages_message" VALUES(167,2,163,0,1,'2010-02-26 03:33:37.903000');
INSERT INTO "messages_message" VALUES(168,23,164,0,1,'2010-02-26 03:41:22.310000');
INSERT INTO "messages_message" VALUES(169,2,164,0,1,'2010-02-26 03:41:22.310000');
INSERT INTO "messages_message" VALUES(170,25,164,0,1,'2010-02-26 03:41:22.310000');
INSERT INTO "messages_message" VALUES(171,23,165,0,1,'2010-02-26 03:50:25.157000');
INSERT INTO "messages_message" VALUES(172,2,165,0,1,'2010-02-26 03:50:25.157000');
INSERT INTO "messages_message" VALUES(173,2,166,0,1,'2010-02-26 03:52:10.690000');
INSERT INTO "messages_message" VALUES(174,10,166,0,1,'2010-02-26 03:52:10.690000');
INSERT INTO "messages_message" VALUES(175,14,166,0,1,'2010-02-26 03:52:10.690000');
INSERT INTO "messages_message" VALUES(176,23,166,0,1,'2010-02-26 03:52:10.690000');
INSERT INTO "messages_message" VALUES(177,23,167,0,1,'2010-02-26 15:51:33.997000');
INSERT INTO "messages_message" VALUES(178,2,167,0,1,'2010-02-26 15:51:33.997000');
INSERT INTO "messages_message" VALUES(179,33,167,0,1,'2010-02-26 15:51:33.997000');
INSERT INTO "messages_message" VALUES(180,23,168,0,1,'2010-02-26 16:07:14.700000');
INSERT INTO "messages_message" VALUES(181,31,168,0,1,'2010-02-26 16:07:14.700000');
INSERT INTO "messages_message" VALUES(182,2,168,0,1,'2010-02-26 16:07:14.700000');
INSERT INTO "messages_message" VALUES(183,23,170,0,1,'2010-02-26 19:15:16.893000');
INSERT INTO "messages_message" VALUES(184,31,170,0,1,'2010-02-26 19:15:16.893000');
INSERT INTO "messages_message" VALUES(185,2,170,0,1,'2010-02-26 19:15:16.893000');
INSERT INTO "messages_message" VALUES(186,38,170,0,1,'2010-02-26 19:15:16.893000');
INSERT INTO "messages_message" VALUES(187,30,171,0,1,'2010-02-26 19:47:41.697000');
INSERT INTO "messages_message" VALUES(188,30,172,0,1,'2010-02-26 19:54:03.447000');
INSERT INTO "messages_message" VALUES(189,42,172,0,1,'2010-02-26 19:54:03.447000');
INSERT INTO "messages_message" VALUES(190,30,173,0,1,'2010-02-26 19:56:22.753000');
INSERT INTO "messages_message" VALUES(191,42,173,0,1,'2010-02-26 19:56:22.753000');
INSERT INTO "messages_message" VALUES(192,23,173,0,1,'2010-02-26 19:56:22.753000');
INSERT INTO "messages_message" VALUES(193,30,174,0,1,'2010-02-26 20:26:41.150000');
INSERT INTO "messages_message" VALUES(194,23,174,0,1,'2010-02-26 20:26:41.150000');
INSERT INTO "messages_message" VALUES(195,2,174,0,1,'2010-02-26 20:26:41.150000');
INSERT INTO "messages_message" VALUES(196,2,175,0,1,'2010-02-26 20:44:01.490000');
INSERT INTO "messages_message" VALUES(197,10,175,0,1,'2010-02-26 20:44:01.490000');
INSERT INTO "messages_message" VALUES(198,14,175,0,1,'2010-02-26 20:44:01.490000');
INSERT INTO "messages_message" VALUES(199,23,175,0,1,'2010-02-26 20:44:01.490000');
INSERT INTO "messages_message" VALUES(200,33,175,0,1,'2010-02-26 20:44:01.490000');
INSERT INTO "messages_message" VALUES(201,44,178,0,1,'2010-02-26 22:15:18.997000');
INSERT INTO "messages_message" VALUES(202,23,180,0,1,'2010-02-26 22:58:09.550000');
INSERT INTO "messages_message" VALUES(203,2,181,0,1,'2010-02-26 23:03:52.940000');
INSERT INTO "messages_message" VALUES(204,23,182,0,1,'2010-02-26 23:04:53.577000');
INSERT INTO "messages_message" VALUES(205,2,182,0,1,'2010-02-26 23:04:53.577000');
INSERT INTO "messages_message" VALUES(206,23,183,0,1,'2010-02-26 23:16:41.967000');
INSERT INTO "messages_message" VALUES(207,2,183,0,1,'2010-02-26 23:16:41.967000');
INSERT INTO "messages_message" VALUES(208,26,183,0,1,'2010-02-26 23:16:41.967000');
INSERT INTO "messages_message" VALUES(209,23,184,0,1,'2010-02-27 02:16:31.970000');
INSERT INTO "messages_message" VALUES(210,2,184,0,1,'2010-02-27 02:16:31.970000');
INSERT INTO "messages_message" VALUES(211,47,184,0,1,'2010-02-27 02:16:31.970000');
INSERT INTO "messages_message" VALUES(212,23,185,0,1,'2010-02-27 02:27:32.303000');
INSERT INTO "messages_message" VALUES(213,31,185,0,1,'2010-02-27 02:27:32.303000');
INSERT INTO "messages_message" VALUES(214,2,185,0,1,'2010-02-27 02:27:32.303000');
INSERT INTO "messages_message" VALUES(215,38,185,0,1,'2010-02-27 02:27:32.303000');
INSERT INTO "messages_message" VALUES(216,40,185,0,1,'2010-02-27 02:27:32.303000');
INSERT INTO "messages_message" VALUES(217,23,186,0,1,'2010-02-27 02:44:04.547000');
INSERT INTO "messages_message" VALUES(218,31,186,0,1,'2010-02-27 02:44:04.547000');
INSERT INTO "messages_message" VALUES(219,2,186,0,1,'2010-02-27 02:44:04.547000');
INSERT INTO "messages_message" VALUES(220,38,186,0,1,'2010-02-27 02:44:04.547000');
INSERT INTO "messages_message" VALUES(221,40,186,0,1,'2010-02-27 02:44:04.547000');
INSERT INTO "messages_message" VALUES(222,24,186,0,1,'2010-02-27 02:44:04.547000');
INSERT INTO "messages_message" VALUES(223,30,187,0,1,'2010-02-28 03:00:56.140000');
INSERT INTO "messages_message" VALUES(224,42,187,0,1,'2010-02-28 03:00:56.140000');
INSERT INTO "messages_message" VALUES(225,23,187,0,1,'2010-02-28 03:00:56.140000');
INSERT INTO "messages_message" VALUES(226,2,187,0,1,'2010-02-28 03:00:56.140000');
INSERT INTO "messages_message" VALUES(227,23,189,0,1,'2010-03-01 20:06:55.973000');
INSERT INTO "messages_message" VALUES(228,23,191,0,1,'2010-03-01 20:41:28.687000');
INSERT INTO "messages_message" VALUES(229,2,191,0,1,'2010-03-01 20:41:28.687000');
INSERT INTO "messages_message" VALUES(230,2,193,0,1,'2010-03-01 20:51:52.190000');
INSERT INTO "messages_message" VALUES(231,5,193,0,1,'2010-03-01 20:51:52.190000');
INSERT INTO "messages_message" VALUES(232,10,194,0,1,'2010-03-01 22:07:15.283000');
INSERT INTO "messages_message" VALUES(233,23,195,0,1,'2010-03-01 22:32:26.107000');
INSERT INTO "messages_message" VALUES(234,2,195,0,1,'2010-03-01 22:32:26.107000');
INSERT INTO "messages_message" VALUES(235,47,195,0,1,'2010-03-01 22:32:26.107000');
INSERT INTO "messages_message" VALUES(236,53,195,0,1,'2010-03-01 22:32:26.107000');
INSERT INTO "messages_message" VALUES(237,10,196,0,1,'2010-03-01 22:52:22.510000');
INSERT INTO "messages_message" VALUES(238,23,196,0,1,'2010-03-01 22:52:22.510000');
INSERT INTO "messages_message" VALUES(239,16,197,0,1,'2010-03-01 23:05:13.920000');
INSERT INTO "messages_message" VALUES(240,2,197,0,1,'2010-03-01 23:05:13.920000');
INSERT INTO "messages_message" VALUES(241,19,197,0,1,'2010-03-01 23:05:13.920000');
INSERT INTO "messages_message" VALUES(242,23,197,0,1,'2010-03-01 23:05:13.920000');
INSERT INTO "messages_message" VALUES(243,16,198,0,1,'2010-03-02 01:51:25.980000');
INSERT INTO "messages_message" VALUES(244,2,198,0,1,'2010-03-02 01:51:25.980000');
INSERT INTO "messages_message" VALUES(245,19,198,0,1,'2010-03-02 01:51:25.980000');
INSERT INTO "messages_message" VALUES(246,23,198,0,1,'2010-03-02 01:51:25.980000');
INSERT INTO "messages_message" VALUES(247,59,198,0,1,'2010-03-02 01:51:25.980000');
INSERT INTO "messages_message" VALUES(248,10,200,0,1,'2010-03-02 02:06:18.193000');
INSERT INTO "messages_message" VALUES(249,10,200,0,1,'2010-03-02 02:06:18.193000');
INSERT INTO "messages_message" VALUES(250,2,201,0,1,'2014-05-09 14:39:20.047000');
INSERT INTO "messages_message" VALUES(251,101,202,0,1,'2014-05-09 14:39:23.890000');
CREATE TABLE "django_admin_log" (
    "id" integer NOT NULL PRIMARY KEY,
    "action_time" datetime NOT NULL,
    "user_id" integer NOT NULL,
    "content_type_id" integer REFERENCES "django_content_type" ("id"),
    "object_id" text,
    "object_repr" varchar(200) NOT NULL,
    "action_flag" smallint unsigned NOT NULL,
    "change_message" text NOT NULL
);
CREATE TABLE "django_flatpage_sites" (
    "id" integer NOT NULL PRIMARY KEY,
    "flatpage_id" integer NOT NULL,
    "site_id" integer NOT NULL REFERENCES "django_site" ("id"),
    UNIQUE ("flatpage_id", "site_id")
);
INSERT INTO "django_flatpage_sites" VALUES(1,1,1);
INSERT INTO "django_flatpage_sites" VALUES(2,2,1);
INSERT INTO "django_flatpage_sites" VALUES(3,3,1);
INSERT INTO "django_flatpage_sites" VALUES(4,4,1);
INSERT INTO "django_flatpage_sites" VALUES(5,5,1);
CREATE TABLE "django_flatpage" (
    "id" integer NOT NULL PRIMARY KEY,
    "url" varchar(100) NOT NULL,
    "title" varchar(200) NOT NULL,
    "content" text NOT NULL,
    "enable_comments" bool NOT NULL,
    "template_name" varchar(70) NOT NULL,
    "registration_required" bool NOT NULL
);
INSERT INTO "django_flatpage" VALUES(1,'/info/faq/','Faq','<h2>Frequently Asked Questions</h2><h3>Contact</h3><p>Contact email: <a href="mailto:admin@biostars.org">admin@biostars.org</a></p><h3>Moderation guidelines</h3><p>Users posting content that does not belong to the site will be notified and required to edit their content. Users may    post commercially motivated posts to the Forum section as long as the topic aligns with the main focus of this site.    Users posting obvious spam will be immediately suspended.</p><h3>User reputation</h3><p>The number next to a user&#39;s name is the sum of upvotes and accepted answers that user has collected.</p><h3>Becoming a moderator</h3><p>Active users above a certain reputation threshold are periodically promoted to moderators. You may also explicitly    ask    for moderation rights or suggest good candidates on the newsgroup. Inactive users that do not visit the site for    extended periods of time may lose their moderation rights.</p><h3>Support for Biostar</h3><p>    Biostar has been developed as an open source software and    has been released with the <b>MIT licence</b>    thanks to the support from the following institutions:<ul>    <li><a href="http://www.psu.edu/">The Pennsylvania State University</a></li>    <li><a href="http://www.nih.gov/">National Institutes of Health (NIH)</a>, specifially grant&nbsp;<a            href="http://ged.msu.edu/downloads/2010-ngs-course-nih-r25.pdf">NIH 5R25HG006243-02, Analyzing Next        Generation        Sequencing Data</a>&nbsp; principal investigator: &nbsp;<a href="http://ged.msu.edu/">Titus Brown (Michigan        State University)</a>,&nbsp;</li></ul></p>',0,'',0);
INSERT INTO "django_flatpage" VALUES(2,'/info/about/','About','<h3>About Biostar</h3>

<p>This site&#39;s focus is&nbsp;<strong>bioinformatics</strong>,&nbsp;<strong>computational genomics</strong>&nbsp;and&nbsp;<strong>biological
    data analysis</strong>. We welcome posts that are:</p>

<ul>
    <li>detailed and specific, written clearly and simply</li>
</ul>

<p>No question is too trivial or too &quot;newbie&quot;.</p>
<p>But we recommend that you make use of the search services to see if your question has already been asked (perhaps
    even answered!) before you ask. But if you end up asking a question that has been asked before, that is fine too.
    Other users will hopefully edit in links to related or similar questions to help future visitors find their way.</p>

<h3>Contact</h3>

<p>To contact the site managers please email&nbsp;<a href="mailto:admin@biostars.org">admin@biostars.org</a>.</p>

<h3>Licensing</h3>

<p>
    <div class="pull-right" style="margin-left:10px"><a rel="license" href="http://creativecommons.org/licenses/by/4.0/">
        <img alt="Creative Commons License"
             style="border-width:0"
             src="http://i.creativecommons.org/l/by/4.0/88x31.png"/></a>
    </div>
    All content on Biostar is licensed via the
    <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons
        Attribution 4.0 International
        License</a>. This means that you should attribute the work
    either to the author and/or to the site, depending on the scope and presentation of the information.
    In addition our community supports the <a href="http://en.wikipedia.org/wiki/Fair_use">fair use policy</a>
    when it comes to content created by users of this site.
</p>


<h3>Citing Biostar</h3>

<p>Parnell LD, Lindenbaum P, Shameer K, Dall&#39;Olio GM, Swan DC, et al. 2011&nbsp;<strong>BioStar: An Online Question
    &amp; Answer Resource for the Bioinformatics Community.</strong>&nbsp;<em>PLoS Comput Biol 7(10)</em>:&nbsp;<a
href="http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002216">e1002216.
    doi:10.1371/journal.pcbi.1002216</a></p>



<h3>Source code</h3>

Biostar source code at <a href="https://github.com/ialbert/biostar-central">biostar-central repository</a>.
Report bugs or feature requests in the <a href="https://github.com/ialbert/biostar-central/issues"> issue tracker</a>.
Copyright by the <a href="https://github.com/ialbert/biostar-central/contributors">BioStar team</a>.',0,'',0);
INSERT INTO "django_flatpage" VALUES(3,'/info/help/','Help','<h2>Help</h2>

<p>&nbsp;</p>

<p>How to get help around here!</p>
',0,'',0);
INSERT INTO "django_flatpage" VALUES(4,'/info/policy/','Policy','<h2>Policy</h2>

<h3>User agreement</h3>

<p>Be well, do good work and keep in touch.</p>

<h3>Privay policy</h3>

<p>We do not give out or sell your information (including name and email address) to anyone.&nbsp;</p>

<h3>Cookies</h3>

<p>The site will set session cookies whenever you visit the site.
    If you do not intend to ever log in, you may deny this
    cookie, but you cannot log in without it.&nbsp;
    <span style="line-height:1.6">Cookies may be also set when you log in, to be able to
        display new posts and messages that you may have received.
        Cookies may also be set by Google or other companies that help us aggregate
        visitor statistics, traffic and other information.</span>
</p>

<h3>Contact</h3>

Contact email: <a href="mailto:admin@biostars.org">admin@biostars.org</a>
',0,'',0);
INSERT INTO "django_flatpage" VALUES(5,'/info/api/','Api','<h2>Biostar API</h2>

<p>
    This is the documentation for Biostar API. If you have additional questions, or believe you have
    encountered a bug, don''t hesitate to post a question on Biostar.
</p>

<h3>General</h3>

<p>
    All API responses are <a href="http://en.wikipedia.org/wiki/JSON">JSON</a>.<br />
    <br />
    Some API responses are cached. Polling for changes should be done sparingly in any case, and
    polling at a rate faster than once a minute (for semantically identical requests) is considered
    abusive.<br />
    <br />
    A number of methods in the Biostar API accept dates as parameters and return dates as properties,
    the format of these dates is documented above. As a general rule, full dates use
    <a href="http://en.wikipedia.org/wiki/ISO_8601">ISO 8601</a> and timestamps are
    in <a href="http://en.wikipedia.org/wiki/Unix_time">unix epoch time</a>.
</p>

<h3>Methods</h3>

<h4>Traffic</h4>
<pre>GET /api/traffic/</pre>
<p>
    Number of post views over the last 60 min filtered by unique IPs.<br />
    <br />
    <b>Fields in response</b><br />
    <ul>
        <li><em>date</em>: the current date, ISO 8601 format.</li>
        <li><em>post_views_last_60_min</em>: number of post views over the last 60 min filtered
        by unique IPs.</li>
        <li><em>timestamp</em>: the current date, unix epoch time format.</li>
    </ul>
    <b>Example</b><br />
    <a href="/api/traffic/">/api/traffic/</a><br />
    <pre>
{
    "date": "2014-05-29T14:59:55.788069",
    "post_views_last_60_min": 850,
    "timestamp": 1401375595
}
    </pre>
</p>

<h4>User</h4>
<pre>GET /api/user/{id}/</pre>
<p>
    General info about a user.<br />
    <br />
    <b>Parameters</b><br />
    <ul>
        <li>
            <em>id</em>: the identifier of the user, a number.
        </li>
    </ul>
    <b>Fields in response</b><br />
    <ul>
        <li><em>date_joined</em>: the date the user joined the website, ISO 8601 format.</li>
        <li><em>id</em>: the identifier of the user, a number.</li>
        <li><em>joined_days_ago</em>: the date the user joined the website, as the number of days
        ago.</li>
        <li><em>last_login</em>: the date of the last login of the user, ISO 8601 format.</li>
        <li><em>name</em>: the name of the user.</li>
        <li><em>vote_count</em>: the number of votes given by the user.</li>
    </ul>
    <b>Example</b><br />
    <a href="/api/user/23/">/api/user/23/</a><br />
    <pre>
{
    "date_joined": "2010-01-18T21:43:55.253000+00:00",
    "id": 23,
    "joined_days_ago": 1614,
    "last_login": "2011-11-08T19:37:21.753000+00:00",
    "name": "Giovanni M Dall''Olio",
    "vote_count": 37
}
    </pre>
</p>

<h4>Post</h4>
<pre>GET /api/post/{id}/</pre>
<p>
    General info about a post.<br />
    <br />
    <b>Parameters</b><br />
    <ul>
        <li>
            <em>id</em>: the identifier of the post, a number.
        </li>
    </ul>
    <b>Fields in response</b><br />
    <ul>
        <li><em>answer_count</em>: number of answers.</li>
        <li><em>author</em>: author name.</li>
        <li><em>author_id</em>: author''s identifier, a number.</li>
        <li><em>book_count</em>: number of bookmarks.</li>
        <li><em>comment_count</em>: number of comments.</li>
        <li><em>creation_date</em>: creation date, ISO 8601 format.</li>
        <li><em>has_accepted</em>: true if the question has an accepted answer, boolean.</li>
        <li><em>id</em>: identifier of the post, a number.</li>
        <li><em>lastedit_date</em>: date of last edit, ISO 8601 format.</li>
        <li><em>lastedit_user_id</em>: user who last edited this post.</li>
        <li><em>parent_id</em>: identifier of the parent post.</li>
        <li><em>rank</em>: rank, a number.</li>
        <li><em>reply_count</em>: number of replies.</li>
        <li><em>root_id</em>: identifier of the root post.</li>
        <li><em>status</em>: status message.</li>
        <li><em>status_id</em>: status'' identifier, a number.</li>
        <li><em>subs_count</em>: number of subscribers following this post.</li>
        <li><em>tag_val</em>: tags.</li>
        <li><em>thread_score</em>: thread''s score.</li>
        <li><em>title</em>: title.</li>
        <li><em>type</em>: type of post.</li>
        <li><em>type_id</em>: type''s identifier for this post.</li>
        <li><em>url</em>: url.</li>
        <li><em>view_count</em>: number of views.</li>
        <li><em>vote_count</em>: number of votes.</li>
        <li><em>xhtml</em>: content.</li>
    </ul>
    <b>Example</b><br />
    <a href="/api/post/25/">/api/post/25/</a><br />
    <pre>
{
    "answer_count": 2,
    "author": "Gue Su",
    "author_id": 18,
    "book_count": 0,
    "comment_count": 0,
    "creation_date": "2009-12-01T20:57:35.300000+00:00",
    "has_accepted": false,
    "id": 25,
    "lastedit_date": "2009-12-01T20:57:35.300000+00:00",
    "lastedit_user_id": 18,
    "parent_id": 24,
    "rank": 0.0,
    "reply_count": 0,
    "root_id": 24,
    "status": "Open",
    "status_id": 1,
    "subs_count": 0,
    "tag_val": "",
    "thread_score": 0,
    "title": "A: How To Set Shrimp Parameters For Best Sensitivity With 35Bp Colorspace Data?",
    "type": "Answer",
    "type_id": 1,
    "url": "http://localhost:8080/p/24/#25",
    "view_count": 0,
    "vote_count": 2,
    "xhtml": "<p>I just read the SHRiMP manual again, but I think that their explanation about -M option may not be enough to answer your question. I usually use the \"seed\" mode by using -s, -n, and -w and the option -M is a new feature of the version 1.3.1, which I have never tried before.</p>\n\n<p>I recommend for you to use the \"seed\" mode--the default would be good, but please adjust the -s option if you want more sensitivity. Always fast speed compensates sensitivity and the -M option seems to exist for this purpose.</p>\n\n<p>Hope my message to be helpful for your project.</p>\n"
}
    </pre>
</p>

<h4>Vote</h4>
<pre>GET /api/vote/{id}/</pre>
<p>
    General info about a vote.<br />
    <br />
    <b>Parameters</b><br />
    <ul>
        <li>
            <em>id</em>: the identifier of the vote, a number.
        </li>
    </ul>
    <b>Fields in response</b><br />
    <ul>
        <li><em>author</em>: author name.</li>
        <li><em>author_id</em>: author''s identifier, a number.</li>
        <li><em>date</em>: date of the vote, ISO 8601 format.</li>
        <li><em>id</em>: identifier of the vote, a number.</li>
        <li><em>post_id</em>: identifier of the voted post.</li>
        <li><em>type</em>: type of vote.</li>
        <li><em>type_id</em>: type''s identifier for this vote.</li>
    </ul>
    <b>Example</b><br />
    <a href="/api/vote/21/">/api/vote/21/</a><br />
    <pre>
{
    "author": "Zhaorong",
    "author_id": 14,
    "date": "2014-04-29T15:02:17.740000+00:00",
    "id": 21,
    "post_id": 26,
    "type": "Upvote",
    "type_id": 0
}
    </pre>
</p>

<h4>Statistics on the <i>N</i>th day</h4>
<pre>GET /api/stats/day/{day}/</pre>
<p>
    Statistics as of the <i>N</i>th day after day-0 (the day of the first ever post).<br />
    <br />
    <b>Parameters</b><br />
    <ul>
        <li>
            <em>day</em>: number of days after day-0, a number.
        </li>
    </ul>
    <b>Fields in response</b><br />
    <ul>
        <li><em>answers</em>: total number of answers as of the given day.</li>
        <li><em>comments</em>: total number of comments as of the given day.</li>
        <li><em>date</em>: date, ISO 8601 format.</li>
        <li><em>new_posts</em>: number of new posts in the given day.</li>
        <li><em>new_users</em>: number of new users in the given day.</li>
        <li><em>new_votes</em>: number of new votes in the given day.</li>
        <li><em>questions</em>: total number of questions as of the given day.</li>
        <li><em>timestamp</em>: date, unix epoch time format.</li>
        <li><em>toplevel</em>: total number of toplevel post as of the given day.</li>
        <li><em>users</em>: total number of users as of the given day.</li>
        <li><em>votes</em>: total number of votes as of the given day.</li>
    </ul>
    <b>Example</b><br />
    <a href="/api/stats/day/5/">/api/stats/day/5/</a><br />
    <pre>
{
    "answers": 6,
    "comments": 0,
    "date": "2009-10-05T00:00:00",
    "new_posts": [
        10,
        11,
        12
    ],
    "new_users": [
        10,
        11
    ],
    "new_votes": [],
    "questions": 6,
    "timestamp": 1254700800,
    "toplevel": 6,
    "users": 10,
    "votes": 0
}
    </pre>
</p>

<h4>Statistics on a date</h4>
<pre>GET /api/stats/date/{year}/{month}/{day}/</pre>
<p>
    Statistics as of the given date.<br />
    <br />
    <b>Parameters</b><br />
    <ul>
        <li><em>year</em>: a number, 4 digits.</li>
        <li><em>month</em>: a number, 2 digits.</li>
        <li><em>day</em>: a number, 2 digits.</li>
    </ul>
    <b>Fields in response</b><br />
    <ul>
        <li><em>answers</em>: total number of answers as of the given date.</li>
        <li><em>comments</em>: total number of comments as of the given date.</li>
        <li><em>date</em>: date, ISO 8601 format.</li>
        <li><em>new_posts</em>: number of new posts in the given date.</li>
        <li><em>new_users</em>: number of new users in the given date.</li>
        <li><em>new_votes</em>: number of new votes in the given date.</li>
        <li><em>questions</em>: total number of questions as of the given date.</li>
        <li><em>timestamp</em>: date, unix epoch time format.</li>
        <li><em>toplevel</em>: total number of toplevel post as of the given date.</li>
        <li><em>users</em>: total number of users as of the given date.</li>
        <li><em>votes</em>: total number of votes as of the given date.</li>
    </ul>
    <b>Example</b><br />
    <a href="/api/stats/date/2009/10/06/">/api/stats/date/2009/10/06/</a><br />
    <pre>
{
    "answers": 9,
    "comments": 0,
    "date": "2009-10-06T00:00:00",
    "new_posts": [
        13,
        14,
        15,
        16
    ],
    "new_users": [
        12,
        13
    ],
    "new_votes": [],
    "questions": 7,
    "timestamp": 1254787200,
    "toplevel": 7,
    "users": 12,
    "votes": 0
}
    </pre>
</p>',0,'',0);
CREATE TABLE "django_session" (
    "session_key" varchar(40) NOT NULL PRIMARY KEY,
    "session_data" text NOT NULL,
    "expire_date" datetime NOT NULL
);
CREATE TABLE "account_emailaddress" (
    "id" integer NOT NULL PRIMARY KEY,
    "user_id" integer NOT NULL,
    "email" varchar(75) NOT NULL UNIQUE,
    "verified" bool NOT NULL,
    "primary" bool NOT NULL
);
CREATE TABLE "account_emailconfirmation" (
    "id" integer NOT NULL PRIMARY KEY,
    "email_address_id" integer NOT NULL REFERENCES "account_emailaddress" ("id"),
    "created" datetime NOT NULL,
    "sent" datetime,
    "key" varchar(64) NOT NULL UNIQUE
);
CREATE TABLE "south_migrationhistory" (
    "id" integer NOT NULL PRIMARY KEY,
    "app_name" varchar(255) NOT NULL,
    "migration" varchar(255) NOT NULL,
    "applied" datetime NOT NULL
);
INSERT INTO "south_migrationhistory" VALUES(1,'users','0001_initial','2015-01-22 16:14:11.363467');
INSERT INTO "south_migrationhistory" VALUES(2,'users','0002_auto__del_field_user_full_score__add_field_user_activity','2015-01-22 16:14:11.389855');
INSERT INTO "south_migrationhistory" VALUES(3,'users','0003_auto__add_tag__add_field_profile_twitter_id__add_field_profile_watch_t','2015-01-22 16:14:11.422613');
INSERT INTO "south_migrationhistory" VALUES(4,'users','0004_auto__add_field_profile_daily_digest__add_field_profile_weekly_digest','2015-01-22 16:14:11.444862');
INSERT INTO "south_migrationhistory" VALUES(5,'users','0005_add_weekly_digest','2015-01-22 16:14:11.453342');
INSERT INTO "south_migrationhistory" VALUES(6,'users','0006_auto__add_field_profile_opt_in','2015-01-22 16:14:11.469735');
INSERT INTO "south_migrationhistory" VALUES(7,'users','0007_auto__del_field_profile_weekly_digest__del_field_profile_daily_digest_','2015-01-22 16:14:11.496959');
INSERT INTO "south_migrationhistory" VALUES(8,'posts','0001_initial','2015-01-22 16:14:12.045255');
INSERT INTO "south_migrationhistory" VALUES(9,'posts','0002_auto__add_data','2015-01-22 16:14:12.114006');
INSERT INTO "south_migrationhistory" VALUES(10,'posts','0003_auto__add_foo','2015-01-22 16:14:12.139514');
INSERT INTO "south_migrationhistory" VALUES(11,'posts','0004_auto__del_data__del_foo__add_emailentry__add_emailsub','2015-01-22 16:14:12.432345');
INSERT INTO "south_migrationhistory" VALUES(12,'badges','0001_initial','2015-01-22 16:14:12.955090');
INSERT INTO "south_migrationhistory" VALUES(13,'badges','0002_auto__del_field_badge_secret__del_field_badge_description__add_field_b','2015-01-22 16:14:12.983829');
INSERT INTO "south_migrationhistory" VALUES(14,'badges','0003_auto__add_field_award_context','2015-01-22 16:14:12.998302');
INSERT INTO "south_migrationhistory" VALUES(15,'planet','0001_initial','2015-01-22 16:14:13.040293');
INSERT INTO "south_migrationhistory" VALUES(16,'planet','0002_auto__add_field_blog_list_order','2015-01-22 16:14:13.051325');
INSERT INTO "south_migrationhistory" VALUES(17,'socialaccount','0001_initial','2015-01-22 16:14:13.089164');
INSERT INTO "south_migrationhistory" VALUES(18,'socialaccount','0002_genericmodels','2015-01-22 16:14:13.147037');
INSERT INTO "south_migrationhistory" VALUES(19,'socialaccount','0003_auto__add_unique_socialaccount_uid_provider','2015-01-22 16:14:13.159355');
INSERT INTO "south_migrationhistory" VALUES(20,'socialaccount','0004_add_sites','2015-01-22 16:14:13.179268');
INSERT INTO "south_migrationhistory" VALUES(21,'socialaccount','0005_set_sites','2015-01-22 16:14:13.192266');
INSERT INTO "south_migrationhistory" VALUES(22,'socialaccount','0006_auto__del_field_socialapp_site','2015-01-22 16:14:13.212295');
INSERT INTO "south_migrationhistory" VALUES(23,'socialaccount','0007_auto__add_field_socialapp_client_id','2015-01-22 16:14:13.243548');
INSERT INTO "south_migrationhistory" VALUES(24,'socialaccount','0008_client_id','2015-01-22 16:14:13.256143');
INSERT INTO "south_migrationhistory" VALUES(25,'socialaccount','0009_auto__add_field_socialtoken_expires_at','2015-01-22 16:14:13.278004');
INSERT INTO "south_migrationhistory" VALUES(26,'socialaccount','0010_auto__chg_field_socialtoken_token','2015-01-22 16:14:13.298224');
INSERT INTO "south_migrationhistory" VALUES(27,'socialaccount','0011_auto__chg_field_socialtoken_token','2015-01-22 16:14:13.319564');
INSERT INTO "south_migrationhistory" VALUES(28,'socialaccount','0012_auto__chg_field_socialtoken_token_secret','2015-01-22 16:14:13.339302');
INSERT INTO "south_migrationhistory" VALUES(29,'djcelery','0001_initial','2015-01-22 16:14:13.421454');
INSERT INTO "south_migrationhistory" VALUES(30,'djcelery','0002_v25_changes','2015-01-22 16:14:13.456191');
INSERT INTO "south_migrationhistory" VALUES(31,'djcelery','0003_v26_changes','2015-01-22 16:14:13.479234');
INSERT INTO "south_migrationhistory" VALUES(32,'djcelery','0004_v30_changes','2015-01-22 16:14:13.496339');
INSERT INTO "south_migrationhistory" VALUES(33,'django','0001_initial','2015-01-22 16:14:13.548555');
INSERT INTO "south_migrationhistory" VALUES(34,'server','0001_initial','2015-01-22 16:14:13.606715');
CREATE TABLE "users_user" ("status" integer NOT NULL, "is_staff" bool NOT NULL, "name" varchar(255) NOT NULL, "site_id" integer, "is_active" bool NOT NULL, "id" integer PRIMARY KEY, "score" integer NOT NULL, "last_login" datetime NOT NULL, "new_messages" integer NOT NULL, "activity" integer NOT NULL, "is_admin" bool NOT NULL, "password" varchar(128) NOT NULL, "type" integer NOT NULL, "email" varchar(255) NOT NULL UNIQUE, "badges" integer NOT NULL, "flair" varchar(15) NOT NULL);
INSERT INTO "users_user" VALUES(0,1,'Biostar Community',NULL,1,1,0,'2014-05-09 14:39:21.988000',0,0,1,'pbkdf2_sha256$12000$v2BdZDmy8Jfc$yZWegbNUk0PzcZ63NzgXiNg5ihS9q044HUyhQvsf2KY=',2,'0@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,1,'István Albert',NULL,1,2,31,'2014-07-11 14:19:18.313000',0,0,1,'pbkdf2_sha256$12000$9qW2F403MdNC$uy0UGLln94XsxdeCW7p2t6KZjQJr9c++/3GNcmSux70=',2,'1@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Fabio',NULL,1,3,3,'2014-04-29 15:02:12.501000',0,0,0,'pbkdf2_sha256$12000$xqGT7BZNn52L$DF76wdfth9Ilt6EFSFzS6wQBCXIX5Kj1MTIgDACpVAQ=',0,'2@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Jason',NULL,1,4,4,'2014-04-29 15:02:12.511000',0,0,0,'pbkdf2_sha256$12000$Z1gu7bXpsT18$H0PVT9g9wHv1rXF33iI9t1RXLnOfMTuKhrQF87WGT+Q=',0,'3@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Zhenhai Zhang',NULL,1,5,5,'2014-04-29 15:02:12.523000',0,0,0,'pbkdf2_sha256$12000$SwQWheoe4DWR$HvnwzSF9zAbW5aJU6dtrWzNS8J/lFe2zLCpYnWqijeo=',0,'4@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Tom Koerber',NULL,1,6,3,'2014-04-29 15:02:12.533000',0,0,0,'pbkdf2_sha256$12000$IwnuXfYFCYNv$CXiYg6rld/tNvdeV8NqSf0ItJ97GQ49TRqg0VCi/y+s=',0,'5@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Suk211',NULL,1,7,6,'2014-04-29 15:02:12.544000',0,0,0,'pbkdf2_sha256$12000$xaF0RyplEZi4$Whthh+miHAbYU7U2INe2es99IZYVxq7C+jhtPAbLZfQ=',0,'6@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Lemon',NULL,1,8,0,'2014-04-29 15:02:12.557000',0,0,0,'pbkdf2_sha256$12000$SxrhjeItCyQ1$yjqCdyEXdvO51kFixpYqW7TbS7QElOP7nLncTlQdnp4=',0,'7@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Wubin Qu',NULL,1,9,0,'2014-04-29 15:02:12.567000',0,0,0,'pbkdf2_sha256$12000$Vu6vHv0Y8ESx$OVkoCTV63vFwqmtPsCM51/iotQuhTpyXFiZ00/U5ib8=',0,'8@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Question Bot',NULL,1,10,9,'2014-04-29 15:02:12.593000',0,0,0,'pbkdf2_sha256$12000$TYnrt8Pwmqgr$S9uwfgPmu6vu+LlQhaqcOUx0Vj4xP6M2BZuTyEjFwFI=',0,'9@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Reka Albert',NULL,1,11,0,'2014-04-29 15:02:12.604000',0,0,0,'pbkdf2_sha256$12000$DP8rjojyp3ym$brpmrCliVz8cxtK3Q4htFEvp0HiZ5oVCxJNNiTDycQI=',0,'10@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Yang Yang',NULL,1,12,2,'2014-04-29 15:02:12.614000',0,0,0,'pbkdf2_sha256$12000$7LJjneJylACC$jEfYvVuRJkRfxShN4AAg1P1hj5iMVbTZ4dUkqANKQJw=',0,'11@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Gue Su Chang',NULL,1,13,1,'2014-04-29 15:02:12.624000',0,0,0,'pbkdf2_sha256$12000$d3mQ7EHUwv8k$INmswQ2vy3R28ssBIinn57i1+pget6ySCtyl6wBJIhQ=',0,'12@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Zhaorong',NULL,1,14,3,'2014-04-29 15:02:12.636000',0,0,0,'pbkdf2_sha256$12000$UFdLVdfM370I$xRLXWq+VifVfvm1k1WVVstwAsFVGQg9NJ231TrZHHlQ=',0,'13@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Nickey',NULL,1,15,0,'2014-04-29 15:02:12.647000',0,0,0,'pbkdf2_sha256$12000$UGNyzHAtnIIJ$VEihUmGm7q5/t/SJvGH9mQPdJixOSoHfljXmkLPbpME=',0,'14@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Renee',NULL,1,16,4,'2014-04-29 15:02:12.657000',0,0,0,'pbkdf2_sha256$12000$aAZZkAtIP7Jj$0LSKb3p9fj5E7nA9w4LNfZydc/hxVueV9gu7iP5supg=',0,'15@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Yu',NULL,1,17,0,'2014-04-29 15:02:12.667000',0,0,0,'pbkdf2_sha256$12000$3CYqycQ4Q5WR$eUqkXV/SL8aNDpEdQWo3dhP1mVVK6do/l4+I/yCSBVs=',0,'16@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Gue Su',NULL,1,18,2,'2014-04-29 15:02:12.678000',0,0,0,'pbkdf2_sha256$12000$WO7H69aEV4Jf$fnZc5IEgOeJeN2H1YB7FIbXvvUc97uJM6P3Yw1vtH58=',0,'17@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Mohammed Islaih',NULL,1,19,2,'2014-04-29 15:02:12.688000',0,0,0,'pbkdf2_sha256$12000$PacGBxydNG0r$JVZGO6FB/BP0NG7OOEpL5O9SmdYsVLDvaD8q89f9y/0=',0,'18@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Alex Reynolds',NULL,1,20,3,'2014-04-29 15:02:12.699000',0,0,0,'pbkdf2_sha256$12000$XZn3ks7DOKGZ$GDtouW5TBZHUXEKGD16ZhZp4GFbLm2X30iyfOprsluk=',0,'19@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'User 4824',NULL,1,21,0,'2014-04-29 15:02:12.709000',0,0,0,'pbkdf2_sha256$12000$RV0TJKXavHsA$4SMx6RZDGt+IRFM8cfO4iCF0V2Zey6CFs6ZQz5dTZ+8=',0,'20@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'User 8226',NULL,1,22,0,'2014-04-29 15:02:12.720000',0,0,0,'pbkdf2_sha256$12000$TGHwkQ27LkJY$SLTaMKnA7rVy1PlmzwSvcKYgxQdnqAegWLgTBc+3Yb8=',0,'21@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Giovanni M Dall''Olio',NULL,1,23,42,'2014-04-29 15:02:12.730000',0,0,0,'pbkdf2_sha256$12000$jSlUQQpIg4zj$1cdW2BEvDSoQO7VdIixeeh6jseW+FZ4+nts2bL3f+so=',1,'22@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Etal',NULL,1,24,7,'2014-04-29 15:02:12.741000',0,0,0,'pbkdf2_sha256$12000$BcCBGGTRvXC8$3dY0tHjT+vLEQXDwx8zZkVRqLCBJiON1libTle8Z6kg=',0,'23@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Fabio',NULL,1,25,3,'2014-04-29 15:02:12.751000',0,0,0,'pbkdf2_sha256$12000$sU81PW5xrTAc$Cu1cFi8pC7c7JfvKlBG2NR27NZTxVAWGFnQuHcsdhkE=',0,'24@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Nicojo',NULL,1,26,4,'2014-04-29 15:02:12.762000',0,0,0,'pbkdf2_sha256$12000$qDuHhjKVFp2x$o+1KCGLxR1RxvsQLYdD5EOegHqJBJ2IGEFfuEwXWGTU=',0,'25@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Allen Yu',NULL,1,27,2,'2014-04-29 15:02:12.772000',0,0,0,'pbkdf2_sha256$12000$rpHuwl7MGZwL$gjhNqj0Gm8HbYfne58W7JRjL1po8ozsmJIiWlF5EduE=',0,'26@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Curious George',NULL,1,28,1,'2014-04-29 15:02:12.783000',0,0,0,'pbkdf2_sha256$12000$oefum8WgLYIN$1unCCIKuvSLoqZkxQ6em9I/NJHraQ8yLUmi5lhM/vEs=',0,'27@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'User 5034',NULL,1,29,0,'2014-04-29 15:02:12.794000',0,0,0,'pbkdf2_sha256$12000$8midSx7wjeGU$IfZR7wT3MSkPHgacZLi5MWwN81JQv4B6J1p0JndWAn8=',0,'28@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Pierre Lindenbaum',NULL,1,30,14,'2014-04-29 15:02:12.805000',0,0,0,'pbkdf2_sha256$12000$dSeZRnZnMkzj$wv9nFF850xquM59gEVATa/KkwCtK58BRqlGkSxaefuU=',2,'29@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Marcos De Carvalho',NULL,1,31,2,'2014-04-29 15:02:12.818000',0,0,0,'pbkdf2_sha256$12000$rr6pEmRIEIfw$100XyJOPNkv6JQq0s8xtKK8b7SfLwtksX79qIstbvoU=',0,'30@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Gustavo Costa',NULL,1,32,0,'2014-04-29 15:02:12.830000',0,0,0,'pbkdf2_sha256$12000$ucdVBO0w9Q7x$Y/mNxoczSRUjJh4+hoUXAGwWJ2i0SMQ6JHkcI4TWwmg=',0,'31@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'User 6996',NULL,1,33,4,'2014-04-29 15:02:12.843000',0,0,0,'pbkdf2_sha256$12000$gTEz50VxJxaZ$sj9mydOzXDyzOiepgCiwPhCykRX/bP5i4yUutRnrzUY=',0,'32@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Paul J. Davis',NULL,1,34,0,'2014-04-29 15:02:12.854000',0,0,0,'pbkdf2_sha256$12000$cM7kWfGaJeVr$za66sGNevXz7UloqJa/2xypORg0Fm77IOja8jR4xQ/A=',0,'33@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'David Nusinow',NULL,1,35,5,'2014-04-29 15:02:12.865000',0,0,0,'pbkdf2_sha256$12000$rVA0q56nigPV$fsadUORb/D6dv07j4+SIl6X0Wc7yzQUI0E055kOC0Xk=',0,'34@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Brentp',NULL,1,36,0,'2014-04-29 15:02:12.875000',0,0,0,'pbkdf2_sha256$12000$4HL5CgiPM0Ia$UcCiPS0N5ms3GTZcEMG6POoveHPTVvsK4TYDMlUh2UM=',1,'35@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Razor',NULL,1,37,0,'2014-04-29 15:02:12.889000',0,0,0,'pbkdf2_sha256$12000$UiIfoDJ8bguC$FPwE36HFeY8gkRIK7hzWfXq9BEP36r27jO1b8BI171g=',0,'36@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Simon Cockell',NULL,1,38,5,'2014-04-29 15:02:12.900000',0,0,0,'pbkdf2_sha256$12000$XYu644JYsxBc$heTzCxR2nE/JPyMYDGy3fJspZGNjlPI0UOjcfLj+ikQ=',1,'37@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Konrad',NULL,1,39,2,'2014-04-29 15:02:12.913000',0,0,0,'pbkdf2_sha256$12000$0e4xhDyZn6UG$8i4EfdOsZhoYVRSN5ivW7VkANmJ8TU+qkCtnURQsyAE=',0,'38@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Biorelated',NULL,1,40,2,'2014-04-29 15:02:12.924000',0,0,0,'pbkdf2_sha256$12000$UwqlVm0drwGy$ym7Fmfkf2gjJePr3LdGwF0EDlwYGNkboHUcdwI5ZXYU=',0,'39@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Yann Abraham',NULL,1,41,0,'2014-04-29 15:02:12.934000',0,0,0,'pbkdf2_sha256$12000$tS2yyrECUCWR$KkpIVlCmMJsgrXy840Y4NLznKjt+znVNyH6LRtnQd1U=',0,'40@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Fernando Muñiz',NULL,1,42,6,'2014-04-29 15:02:12.948000',0,0,0,'pbkdf2_sha256$12000$UR24VGeEwAYn$rJgvCgqhgKWU7avwI8oqxDrXkSBR3f+dmW61snodUT4=',0,'41@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'User 1402',NULL,1,43,0,'2014-04-29 15:02:12.959000',0,0,0,'pbkdf2_sha256$12000$TAnq0LpvHxen$fBlD33enZF59h6Ow1kh5DZ5T5EvAEGGm44PWvQWtXJY=',0,'42@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Andrew',NULL,1,44,1,'2014-04-29 15:02:12.969000',0,0,0,'pbkdf2_sha256$12000$QSevRLvHj4fs$/EyGRPMnRfQgKJpqvQdG3FrlVNmDN/xyQVs2IqB0V4Q=',0,'43@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Alex',NULL,1,45,1,'2014-04-29 15:02:12.980000',0,0,0,'pbkdf2_sha256$12000$xgfjPCTBuZh4$vwYrG3WezuinPOLeYW7gnTOkxpkXRQTFz0cwzCCvCi8=',0,'44@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Owen',NULL,1,46,0,'2014-04-29 15:02:12.990000',0,0,0,'pbkdf2_sha256$12000$8dLFyZZD7biT$YLgRWfpWkKC9cYoaWWJMbUI5vBCacdIvGsCJDJ4BF7g=',0,'45@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Chris',NULL,1,47,3,'2014-04-29 15:02:13.001000',0,0,0,'pbkdf2_sha256$12000$oMoTXWGRFHP1$hMcQDKHlOrbUu7zFhv6qFxFdKctwFiVxxBm0kt6wsjw=',0,'46@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Kelly O.',NULL,1,48,0,'2014-04-29 15:02:13.013000',0,0,0,'pbkdf2_sha256$12000$9v2mWmQ7WHfC$6gm9/LF+awSRyswDAww4ZscmbMuS9BMN1uQIZgwvYfM=',0,'47@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'User 3704',NULL,1,49,0,'2014-04-29 15:02:13.023000',0,0,0,'pbkdf2_sha256$12000$Kpd2XLAgDrsw$7tf44RV/4+HTOl+93LnTtfXkliIwIEsIXXhX6c6TGkU=',0,'48@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Greg',NULL,1,50,0,'2014-04-29 15:02:13.034000',0,0,0,'pbkdf2_sha256$12000$cfxefsjduUzM$3E/1A3nxt/kKCMFMut+QbKd4V9l1LQs65IQrfRUA+LE=',0,'49@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Pedrobeltrao',NULL,1,51,0,'2014-04-29 15:02:13.046000',0,0,0,'pbkdf2_sha256$12000$SvRMm9waJ3qO$hCkfgGmYn/VIVwvNMc3BnUnMdDFnVa9g76n/u+DEtGA=',0,'50@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Liam Thompson',NULL,1,52,0,'2014-04-29 15:02:13.060000',0,0,0,'pbkdf2_sha256$12000$aOGcLeQEfIPn$sQdkkKJWC72P1zoBmddxLQDxl+3iqUYgHYhXHc0/nCQ=',0,'51@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Michael Barton',NULL,1,53,3,'2014-04-29 15:02:13.072000',0,0,0,'pbkdf2_sha256$12000$FBhI1GpNgZrO$5cu8QS477St6WFya4dIkVurUxHqG2zRSd6PmlHP1Skg=',0,'52@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Schrodinger''S Cat',NULL,1,54,0,'2014-04-29 15:02:13.085000',0,0,0,'pbkdf2_sha256$12000$cHejL2WNI43X$XgstVRHdIe1waNfCkzfanSyMSUm30Gg+Ms+vG1/Pn70=',0,'53@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Michael Dondrup',NULL,1,55,5,'2014-04-29 15:02:13.096000',0,0,0,'pbkdf2_sha256$12000$sZbNI0zC2u60$08s9CquCBfOMyn1w+kc+R4PjxOokxh0E+70xgBo7x3Q=',1,'54@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Brad Chapman',NULL,1,56,0,'2014-04-29 15:02:13.109000',0,0,0,'pbkdf2_sha256$12000$UvODgEajtzmh$TdOr4eJywyT+W4gx5dCze+t4Vll+6/VGEJVvglYD1p0=',1,'55@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Paulati',NULL,1,57,0,'2014-04-29 15:02:13.121000',0,0,0,'pbkdf2_sha256$12000$ZgG779ob86Sl$wl94Aq/66lmyn6YMQgP9AC7HsdjY7pcM9TtPIvT58yQ=',0,'56@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Dave Bridges',NULL,1,58,0,'2014-04-29 15:02:13.134000',0,0,0,'pbkdf2_sha256$12000$ov8ArYAwjJlr$gFUKnAYbkeDj1EvlljR/BN+Zn5kDv+jDHCjJfKSGMCo=',0,'57@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Daniel Swan',NULL,1,59,2,'2014-04-29 15:02:13.144000',0,0,0,'pbkdf2_sha256$12000$tFSFE4jMaOJ8$+PMRrl7qwOac8w1K7q+nJIdHg93uwTBYaMAa9corEJQ=',1,'58@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Lukfor',NULL,1,60,0,'2014-04-29 15:02:13.160000',0,0,0,'pbkdf2_sha256$12000$IXzi1Xn9Czas$1NrrU+G+MTc17HHHbiWfCM+q3C+O+Z4CDzdQnIPQ7os=',0,'59@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Chris Fields',NULL,1,61,0,'2014-04-29 15:02:13.171000',0,0,0,'pbkdf2_sha256$12000$3m6Wp5oeGCz6$rCsr8gIYXKR1TYd32eL+w4tqMair81M8LBimwEqrSs0=',0,'60@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Darked89',NULL,1,62,0,'2014-04-29 15:02:13.184000',0,0,0,'pbkdf2_sha256$12000$x6BrtZZDJNDF$NE75vuPcXt5No5fbfRVCE2WEO1TrmNTpL7ESAODtMow=',1,'61@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Lmartinho',NULL,1,63,0,'2014-04-29 15:02:13.195000',0,0,0,'pbkdf2_sha256$12000$sH5jvrYb2vyn$PZ2A17ssUR0fvm+mq0ThBgB4Ih9EMaXmRSwzCNvynkE=',0,'62@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Manuel Corpas',NULL,1,64,0,'2014-04-29 15:02:13.206000',0,0,0,'pbkdf2_sha256$12000$4UaQ3C3diXy5$SwjuVgnbHwBzmiPFQfSCpct1a7YQIMYXaOuPqaJEwD4=',0,'63@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Abhishek Tiwari',NULL,1,65,0,'2014-04-29 15:02:13.220000',0,0,0,'pbkdf2_sha256$12000$qAFmBJPZgTSh$xgnqaXKc7bPecKLtMUPa/Ea4kXRmiUIaIcB58uuG72c=',0,'64@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Neilfws',NULL,1,66,0,'2014-04-29 15:02:13.245000',0,0,0,'pbkdf2_sha256$12000$3uS5mmrfxsYn$uPZoPv/tgqDF6OCrKeUwxiwfDAQlpRchJwYtvZyiYWA=',2,'65@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Piotr Byzia',NULL,1,67,0,'2014-04-29 15:02:13.267000',0,0,0,'pbkdf2_sha256$12000$EVvOhIOv8eQP$1ygCXtMXg7OzNdvvgKPoF4+jEGAeosXzxRDjNf1E2UE=',0,'66@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Jeroen Van Goey',NULL,1,68,0,'2014-04-29 15:02:13.281000',0,0,0,'pbkdf2_sha256$12000$KD1LVVbEXmTR$FTjcxzC/7zXXG6NLWc1Uh2lN0HWL1SFVIFAue6REWJc=',0,'67@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Vince',NULL,1,69,0,'2014-04-29 15:02:13.296000',0,0,0,'pbkdf2_sha256$12000$QqzZpb87JWKo$AJUgGmiNCwHEBjZxllvCsIU/ZqlZNSkZz733vv5nr9s=',0,'68@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Tom Walsh',NULL,1,70,0,'2014-04-29 15:02:13.306000',0,0,0,'pbkdf2_sha256$12000$6AKEk1eW77tA$ncxfHiq7OJMq9S6ewqtYmt5K6tm6Agt0yLtW53Iawro=',0,'69@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Egon Willighagen',NULL,1,71,0,'2014-04-29 15:02:13.333000',0,0,0,'pbkdf2_sha256$12000$xIbUkiEIEySs$k6+6hPRRD4NEM9xJ0hig0PIf6koGHRtg7QO0DVpQCRw=',1,'70@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Mndoci',NULL,1,72,0,'2014-04-29 15:02:13.349000',0,0,0,'pbkdf2_sha256$12000$ofoMoq3WABGc$cjeS/IWd0EnaEPNejlepL4DqHm3ad0QzRISAffTXEps=',0,'71@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Jeremy Leipzig',NULL,1,73,0,'2014-04-29 15:02:13.365000',0,0,0,'pbkdf2_sha256$12000$5vXXOzn7tPKn$1UQvXBDX6kmkfoypsARwqPYNrPPUB49c6pfpHT+MmKU=',1,'72@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Mmarchin',NULL,1,74,0,'2014-04-29 15:02:13.380000',0,0,0,'pbkdf2_sha256$12000$MGJpwYIMsQOi$I9JWhUxQRmBgcCMP2DdmKKsJtEkEKyQ/3vhsFwyemHQ=',0,'73@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Paolo',NULL,1,75,0,'2014-04-29 15:02:13.391000',0,0,0,'pbkdf2_sha256$12000$92BQz3pwTcOO$d411VgY6cCq7qlfNeJwpsVYZ/JT1kWlASIPkhmDIZa8=',0,'74@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Andrew',NULL,1,76,0,'2014-04-29 15:02:13.404000',0,0,0,'pbkdf2_sha256$12000$gasAz3kX8aXj$ev2WyN9ZaByNLTIOKZEceYNv6i5IaGoUHAcXZHdIjJU=',0,'75@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Eleanor Howe',NULL,1,77,0,'2014-04-29 15:02:13.415000',0,0,0,'pbkdf2_sha256$12000$DUN1fMdrF2SR$5ODQ4g1j7p7UzPEe/6W/CKpUGKTQzRYhfcfIAOMmq34=',0,'76@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Jboveda',NULL,1,78,0,'2014-04-29 15:02:13.434000',0,0,0,'pbkdf2_sha256$12000$FIXFvfZQTS83$kBWJNjwh8cVMkuvqAtG0oo61pknQBCaGMYLQRw1SfP4=',0,'77@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Tim',NULL,1,79,0,'2014-04-29 15:02:13.446000',0,0,0,'pbkdf2_sha256$12000$WveuTv2dX7iJ$YoRcZ18ysiR//FFmii4emIPwyIHwznlZYTtvUt4ngZM=',0,'78@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Daniel Jurczak',NULL,1,80,0,'2014-04-29 15:02:13.458000',0,0,0,'pbkdf2_sha256$12000$NBngO2uzk6y7$wSXv00L+6yWPOzCQ/knWe44Xox+suvhecOWnhiHSWE4=',0,'79@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'User 1073',NULL,1,81,0,'2014-04-29 15:02:13.468000',0,0,0,'pbkdf2_sha256$12000$O8oYZpiqY314$Y6irPwDq8N1/GIgL2lOgLo5tnOUNSVEr7pokJ/clBjM=',0,'80@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Geoffjentry',NULL,1,82,0,'2014-04-29 15:02:13.482000',0,0,0,'pbkdf2_sha256$12000$anjIFGaJSxcl$YLL2Xw+QHLh98mN48ENGLVJUMiyNGj01FAjNnutiing=',0,'81@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'User 3505',NULL,1,83,0,'2014-04-29 15:02:13.495000',0,0,0,'pbkdf2_sha256$12000$2Hwvgl2XRACg$pcBjf8kPDB/ShYNp47laLenjOwTgayb025LOdKWGffE=',0,'82@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Anshu',NULL,1,84,0,'2014-04-29 15:02:13.506000',0,0,0,'pbkdf2_sha256$12000$0gD5k5pLyFMU$VBeYrroiSSoL9+w4n/XGdYMi3PLwTQqyfgrKzT6CLV8=',0,'83@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Bryan Maloney',NULL,1,85,0,'2014-04-29 15:02:13.517000',0,0,0,'pbkdf2_sha256$12000$F49036a6RZu5$WBS83LIKcm535yw2noIhyGKQBqv/ndxk1M9lk8Q2GXk=',0,'84@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Andrew Su',NULL,1,86,0,'2014-04-29 15:02:13.530000',0,0,0,'pbkdf2_sha256$12000$xLzfu0zOdYgl$g/avnMo0wkDs7FKwnSvCMW5Nn0NopGKmN9FPA8rmo+s=',1,'85@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Khader Shameer',NULL,1,87,0,'2014-04-29 15:02:13.542000',0,0,0,'pbkdf2_sha256$12000$f7DHaswt7bGg$dXRcjhjmUCyeZlESq0Hi4Orl+L62CYdZa3Zcnfq5MIM=',1,'86@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Pierre',NULL,1,88,0,'2014-04-29 15:02:13.554000',0,0,0,'pbkdf2_sha256$12000$erfcUO9QwlNF$JRt/ogMgZPkFqf2POlVqqtTFNCzG1PcJ8MBfv+VABKo=',0,'87@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Yuri',NULL,1,89,0,'2014-04-29 15:02:13.575000',0,0,0,'pbkdf2_sha256$12000$cewFu6xQ6Ae1$z95BGoh2E+Yh3Tl3JDWa2+0Tn34XYdvOLvga4ckhJF4=',0,'88@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Allpowerde',NULL,1,90,0,'2014-04-29 15:02:13.588000',0,0,0,'pbkdf2_sha256$12000$7oJ63SMrqkmL$2QilTjsJg7+R8oFbsjFA3TkzLbt6TdvrU4VYe3Dbei0=',0,'89@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Mikael Huss',NULL,1,91,0,'2014-04-29 15:02:13.613000',0,0,0,'pbkdf2_sha256$12000$QjVJyGehvQXY$NCGKqDQirpXrxN67pAng7cw5yMmrpB2BJ2C7NN2Nuao=',0,'90@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Perry',NULL,1,92,0,'2014-04-29 15:02:13.632000',0,0,0,'pbkdf2_sha256$12000$qkRsBi11opyj$Kb3rMqW1yIuA9X6+LxeJtFS5+cdurftyTBGT3jfHvDA=',0,'91@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'User 5106',NULL,1,93,0,'2014-04-29 15:02:13.649000',0,0,0,'pbkdf2_sha256$12000$vTS5D5OcYAny$oztMskBK5ZezMRW9JLWvhuZ1zsLbnuDVtrbi0Tt2WWw=',0,'92@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Nir London',NULL,1,94,0,'2014-04-29 15:02:13.662000',0,0,0,'pbkdf2_sha256$12000$rzaCSA5TUvNf$9+5WN/Egbq8Y8IB4Weoc+yVEKGarBg0PntbJEAT4AXM=',0,'93@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Bioinfo',NULL,1,95,0,'2014-04-29 15:02:13.683000',0,0,0,'pbkdf2_sha256$12000$cqZGG8CKNJQj$qb8AEKQGWLd6p6maQ7zQ8+0QpCEBToYPwsPYJiYid64=',0,'94@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Jc',NULL,1,96,0,'2014-04-29 15:02:13.695000',0,0,0,'pbkdf2_sha256$12000$V4ZiZxoOxTfu$0d2e6vFap5TYkCvpXuLcMINZ4SFXNns2ojxsLPYpXhc=',0,'95@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Michael Hoffman',NULL,1,97,0,'2014-04-29 15:02:13.712000',0,0,0,'pbkdf2_sha256$12000$dNsbB7ejmsyG$E0QUnPMyuHjx48eZjoDiOg2JO+vrdX6TpP3DBWBD+FQ=',0,'96@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Nancy Parmalee',NULL,1,98,0,'2014-04-29 15:02:13.732000',0,0,0,'pbkdf2_sha256$12000$2RaQnkNvUKyH$puGanmv6F/QUs592SbpkRx8yvyQdhXG6PgV3pwGmzkM=',0,'97@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Avilella',NULL,1,99,0,'2014-04-29 15:02:13.743000',0,0,0,'pbkdf2_sha256$12000$2FTVddwp9Hty$ILbwZtEYc50AVQBkQACUYk3d8h3nxmVWl0ufjUQhFQs=',0,'98@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'Ryan',NULL,1,100,0,'2014-04-29 15:02:13.768000',0,0,0,'pbkdf2_sha256$12000$MSN9hQi1fryA$Nhp6n8kOq0HQD40WO5pCWWJyi8IpDdimWgv+ErlLHr0=',0,'99@foo.bar',0,'');
INSERT INTO "users_user" VALUES(0,0,'john',NULL,1,101,0,'2014-05-09 14:39:23.915000',0,0,0,'pbkdf2_sha256$12000$2E3GGcg1ofYT$QeUtTdcst73E0+5l6vtr+ueP1qmL11ca0aKEyJRg6Zc=',0,'100@foo.bar',0,'');
CREATE TABLE "users_emaillist" ("id" integer NOT NULL PRIMARY KEY, "date" datetime NOT NULL, "email" varchar(255) NOT NULL UNIQUE, "type" integer NOT NULL, "active" bool NOT NULL);
CREATE TABLE "users_tag" ("id" integer NOT NULL PRIMARY KEY, "name" text NOT NULL);
INSERT INTO "users_tag" VALUES(1,'galaxy');
INSERT INTO "users_tag" VALUES(2,'minia');
CREATE TABLE "users_profile_tags" ("id" integer NOT NULL PRIMARY KEY, "profile_id" integer NOT NULL, "tag_id" integer NOT NULL);
INSERT INTO "users_profile_tags" VALUES(1,2,1);
INSERT INTO "users_profile_tags" VALUES(2,2,2);
CREATE TABLE "users_profile" ("website" varchar(255) NOT NULL, "info" text, "user_id" integer NOT NULL UNIQUE, "uuid" varchar(255) NOT NULL UNIQUE, "digest_prefs" integer NOT NULL, "opt_in" bool NOT NULL, "twitter_id" varchar(255) NOT NULL, "my_tags" text NOT NULL, "flag" integer NOT NULL, "last_login" datetime NOT NULL, "location" varchar(255) NOT NULL, "watched_tags" varchar(100) NOT NULL, "message_prefs" integer NOT NULL, "scholar" varchar(255) NOT NULL, "id" integer PRIMARY KEY, "date_joined" datetime NOT NULL);
INSERT INTO "users_profile" VALUES('','',1,'9428b06c8c04f25b5abf95d8354b2571',2,0,'','',0,'2014-05-09 14:39:21.999000','State College, USA','',3,'',1,'2014-04-29 15:02:11.528000');
INSERT INTO "users_profile" VALUES('http://www.personal.psu.edu/iua1/','<p>Current position: Research Associate Professor of Bioinformatics at Penn State</p>
',2,'7db961d01bd62e569f584bd0aeb035aa',2,0,'zumm','',0,'2014-07-11 14:19:18.323000','University Park','galaxy, minia',3,'',2,'2009-09-30 19:30:39.570000');
INSERT INTO "users_profile" VALUES('','',3,'2c95a46647f03e33aea79243ec353b26',2,0,'','',0,'2009-10-05 21:43:12.173000','','',3,'',3,'2009-09-30 22:09:06.600000');
INSERT INTO "users_profile" VALUES('','<p>Grad student at PSU</p>
',4,'1ea9fa61da4b44c4329942cdf51710e6',2,0,'','',0,'2011-11-04 20:45:04.227000','','',3,'',4,'2009-09-30 23:10:54.623000');
INSERT INTO "users_profile" VALUES('','',5,'8a92ca71e614f2669e5114cce384bc6a',2,0,'','',0,'2011-07-06 23:44:55.723000','502 Wartik Lab, Penn State Univ','',3,'',5,'2009-10-01 00:48:08.723000');
INSERT INTO "users_profile" VALUES('','',6,'87e4c8d6cfa89ae3bd8616cca02189c4',2,0,'','',0,'2009-10-01 01:32:29.033000','','',3,'',6,'2009-10-01 01:32:29.033000');
INSERT INTO "users_profile" VALUES('http://twitter.com/suk211','<p>interested systems biology,protein folding ,genomics</p>
',7,'b6bc828bfce3913972cf286618113ce9',2,0,'','',0,'2011-07-07 09:32:47.793000','state college','',3,'',7,'2009-10-01 01:35:27.973000');
INSERT INTO "users_profile" VALUES('','',8,'e038b231cc619eef4639e884015333ca',2,0,'','',0,'2010-10-14 05:59:00.327000','China,WH','',3,'',8,'2009-10-03 09:02:07.727000');
INSERT INTO "users_profile" VALUES('http://quwubin.sinaapp.com/','<p>Bioinformatics researcher, author of MFEprimer, MPprimer and CelRNAi.</p>
',9,'8342495beebd4d7c2fb8d3e3f97bf64d',2,0,'','',0,'2011-11-04 06:28:53.817000','China','',3,'',9,'2009-10-03 16:15:52.077000');
INSERT INTO "users_profile" VALUES('','',10,'de65e553cc0085e0162ebc7cbc7204e7',2,0,'','',0,'2013-04-22 18:30:22.763000','','',3,'',10,'2009-10-05 21:44:30.853000');
INSERT INTO "users_profile" VALUES('http://www.phys.psu.edu/~ralbert/','',11,'c09a115704330b6c839291009efb0b64',2,0,'','',0,'2009-10-05 22:47:25.690000','University Park','',3,'',11,'2009-10-05 22:47:25.690000');
INSERT INTO "users_profile" VALUES('','',12,'b1f45396c03358d4c998383fb1510d81',2,0,'','',0,'2009-10-07 07:45:15.080000','','',3,'',12,'2009-10-07 00:58:10.083000');
INSERT INTO "users_profile" VALUES('','',13,'a4700b3df7890aead3bc10749213524d',2,0,'','',0,'2009-10-07 19:04:30.217000','','',3,'',13,'2009-10-07 04:04:41.700000');
INSERT INTO "users_profile" VALUES('','<p>A grad student in bioinformatics studying plant small RNAs.</p>
',14,'397f3990a96766bf70a13b7ac42daaf0',2,0,'','',0,'2011-10-29 04:18:27.637000','State College, PA','',3,'',14,'2009-10-09 08:54:13.480000');
INSERT INTO "users_profile" VALUES('','',15,'b12240f4a599d5f23445e323a142cd5b',2,0,'','',0,'2009-10-18 09:22:53.687000','','',3,'',15,'2009-10-18 09:22:53.687000');
INSERT INTO "users_profile" VALUES('','',16,'fc355f3a153dc09cdf592308ee26fc0f',2,0,'','',0,'2010-02-21 03:13:56.363000','','',3,'',16,'2009-10-23 23:40:27.030000');
INSERT INTO "users_profile" VALUES('','',17,'ff811eb67a6b661c2f26f8d6d5b03975',2,0,'','',0,'2011-10-25 03:24:37.190000','','',3,'',17,'2009-11-05 00:39:12.007000');
INSERT INTO "users_profile" VALUES('','',18,'8920f2b251f1e0f4405777ca86bbf957',2,0,'','',0,'2009-12-01 22:51:54.427000','','',3,'',18,'2009-12-01 20:57:35.020000');
INSERT INTO "users_profile" VALUES('','',19,'b8f90a1927779e2fc1ec2469859420c8',2,0,'','',0,'2009-12-09 03:45:54.220000','','',3,'',19,'2009-12-09 03:45:54.220000');
INSERT INTO "users_profile" VALUES('','',20,'300f2a04cb6f5d8d27e732faad01f12c',2,0,'','',0,'2011-11-08 18:23:58.050000','','',3,'',20,'2009-12-23 07:10:22.660000');
INSERT INTO "users_profile" VALUES('','',21,'2cdc27d640e61b16951f016af7b35aaa',2,0,'','',0,'2009-12-28 19:56:46.033000','','',3,'',21,'2009-12-28 19:56:46.033000');
INSERT INTO "users_profile" VALUES('','',22,'75bd966bfb9d1bf7f254020f169a5863',2,0,'','',0,'2010-04-10 08:15:19.790000','','',3,'',22,'2009-12-31 19:26:29.283000');
INSERT INTO "users_profile" VALUES('http://bioinfoblog.it','<p>PhD student at the Pompeu Fabra University in Barcelona. Interested in bioinformatics, R, python programming and agile techniques.</p>
',23,'412b84490c0ce653cfb8729cc2f9cef8',2,0,'','',0,'2011-11-08 19:37:21.753000','Barcelona, Spain','',3,'',23,'2010-01-18 21:43:55.253000');
INSERT INTO "users_profile" VALUES('http://etalog.blogspot.com/','',24,'c5b6e90e94c47609a7d1bec571083e1d',2,0,'','',0,'2011-11-08 05:26:49.263000','Athens, GA','',3,'',24,'2010-01-23 23:01:47.617000');
INSERT INTO "users_profile" VALUES('','',25,'6b1a47a5753751738f8c739824605036',2,0,'','',0,'2010-01-27 06:38:35.550000','','',3,'',25,'2010-01-27 05:42:13.980000');
INSERT INTO "users_profile" VALUES('','',26,'ca27575af4b79d86d5df1b03ce0dc93f',2,0,'','',0,'2011-10-28 12:58:08.517000','','',3,'',26,'2010-02-13 03:29:24.910000');
INSERT INTO "users_profile" VALUES('','',27,'f86e90b7decb0efc9c119c2cd718c362',2,0,'','',0,'2010-03-11 08:42:46.730000','','',3,'',27,'2010-02-19 13:25:13.603000');
INSERT INTO "users_profile" VALUES('','',28,'7f8f978a6f0b2071fba8cc52ecb4c38f',2,0,'','',0,'2010-02-20 01:35:30.157000','','',3,'',28,'2010-02-19 19:33:45.333000');
INSERT INTO "users_profile" VALUES('','',29,'d24e7138efa7444a7b0ee4fc99200979',2,0,'','',0,'2011-03-17 11:18:47.650000','','',3,'',29,'2010-02-22 09:24:18.933000');
INSERT INTO "users_profile" VALUES('http://plindenbaum.blogspot.com','<p>Virology, Bioinformatics , Genetics</p>

<p><a href="http://twitter.com/yokofakun">http://twitter.com/yokofakun</a></p>

<p>Neil Saunders is an impostor because he doesn''t exist.</p>

<p>Please, don''t vote for him :-)</p>
',30,'15ee06265e393f65a8ec6404df6a9ba9',2,0,'','',0,'2011-11-08 20:35:24.150000','France','',3,'',30,'2010-02-25 23:27:15.560000');
INSERT INTO "users_profile" VALUES('http://www.marcosdecarvalho.com','<p>Doctoral student.</p>
',31,'e872facb6b47883644f0965eb2192f95',2,0,'','',0,'2011-08-12 02:05:42.993000','Porto Alegre, RS, Brasil','',3,'',31,'2010-02-25 23:42:49.150000');
INSERT INTO "users_profile" VALUES('http://www.lge.ibi.unicamp.br','<p>Phd candidate at State University of Campinas (UNICAMP)
Bs in Computer Science (UFS)</p>

<p>Bioinformatician with experience in genomics, transcriptomics, gene finding and next-generation sequencing.</p>
',32,'48ba9d1464512088542c06b1a55a860b',2,0,'','',0,'2010-02-26 00:22:10.283000','Campinas, Brazil','',3,'',32,'2010-02-26 00:22:10.283000');
INSERT INTO "users_profile" VALUES('','',33,'2aa464fbd7c268351271132fa3b13a74',2,0,'','',0,'2010-02-26 06:09:01.270000','','',3,'',33,'2010-02-26 00:28:54.453000');
INSERT INTO "users_profile" VALUES('http://www.davispj.com','',34,'659247b1e4f22d03109bbdf33bacb608',2,0,'','',0,'2011-06-10 02:08:05.940000','Ipswich, MA','',3,'',34,'2010-02-26 00:31:12.473000');
INSERT INTO "users_profile" VALUES('','',35,'b9fb010aac2a4bc44c112cff91312af3',2,0,'','',0,'2011-11-02 04:06:55.930000','Boston, MA','',3,'',35,'2010-02-26 03:30:36.047000');
INSERT INTO "users_profile" VALUES('https://github.com/brentp','<p>I do bioinformatics at a hospital in denver.</p>

<p>
<a href="http://github.com/brentp/" rel="nofollow">some of my code</a>
</p>

<p>
<a href="http://hackmap.blogspot.com/" rel="nofollow">and my blog</a>
</p>

<p>and 
<a href="http://twitter.com/#!/brent_p" rel="nofollow">twitter</a></p>
',36,'ad3bfa7ffb49acbde3f9617eaab24dd0',2,0,'','',0,'2011-11-08 19:32:43.780000','Denver, Colorado','',3,'',36,'2010-02-26 03:33:47.560000');
INSERT INTO "users_profile" VALUES('http://biopunk.hu','',37,'e268c14fe9193ec154cd76924f2a71c3',2,0,'','',0,'2010-06-21 21:27:38.803000','Budapest','',3,'',37,'2010-02-26 13:58:44.137000');
INSERT INTO "users_profile" VALUES('http://fuzzierlogic.com/','<p></p><p>Bioinformatician, support scientist, Python programmer. Full bio on <a href="http://blog.fuzzierlogic.com/about" rel="nofollow">my blog</a>. You can also follow me <a href="http://twitter.com/sjcockell" rel="nofollow">on Twitter</a>.</p>

<p></p><p>I work with <a href="http://biostar.stackexchange.com/users/59/daniel-swan" rel="nofollow">Daniel</a> in the <a href="http://bsu.ncl.ac.uk/" rel="nofollow">Bioinformatics Support Unit</a> of Newcastle University.</p>
',38,'e4a60f08057d7f587542f0b2826472e3',2,0,'','',0,'2011-11-08 16:47:47.327000','Newcastle','',3,'',38,'2010-02-26 15:26:42.033000');
INSERT INTO "users_profile" VALUES('http://konrad.foerstner.org','',39,'959a865854e2d90ebe04bb11440829de',2,0,'','',0,'2011-11-08 20:33:26.513000','','',3,'',39,'2010-02-26 19:01:06.773000');
INSERT INTO "users_profile" VALUES('http://www.biorelated.com','',40,'8ca824b97e40af61ed34dcf60ccdda16',2,0,'','',0,'2011-11-08 01:51:55.693000','','',3,'',40,'2010-02-26 19:06:05.413000');
INSERT INTO "users_profile" VALUES('','<p>computational biology for the pharmaceutical industry, father, geek, climber, biker</p>
',41,'f49922dc6480ab18cdb0b034e0a837d5',2,0,'','',0,'2011-09-05 14:27:23.120000','Basel, Switzerland','',3,'',41,'2010-02-26 19:37:34.387000');
INSERT INTO "users_profile" VALUES('','',42,'ccac10608aa8463c4a841e7bd9a92fb6',2,0,'','',0,'2010-02-26 20:26:41.170000','','',3,'',42,'2010-02-26 19:47:41.540000');
INSERT INTO "users_profile" VALUES('','',43,'f541c1df4528a879a9637c9a1dee57b0',2,0,'','',0,'2010-02-26 20:47:27.600000','','',3,'',43,'2010-02-26 20:47:27.600000');
INSERT INTO "users_profile" VALUES('','',44,'1a1dc99c3e2b299c0c81aa0606e1e8af',2,0,'','',0,'2011-11-07 20:16:10.883000','','',3,'',44,'2010-02-26 21:49:18.470000');
INSERT INTO "users_profile" VALUES('','',45,'0a0e02f753dcf32806420c68471c5fd9',2,0,'','',0,'2010-03-10 04:07:11.900000','','',3,'',45,'2010-02-26 22:06:54.723000');
INSERT INTO "users_profile" VALUES('','',46,'c59f65d8e8305712a786e6c6d1e66cc5',2,0,'','',0,'2011-05-21 06:41:18.377000','','',3,'',46,'2010-02-26 22:54:23.013000');
INSERT INTO "users_profile" VALUES('http://www.rostlab.org/~schaefer/','<p>Researcher in structural bioinformatics</p>
',47,'853ff6827bf2e1b575bae774a234ff4e',2,0,'','',0,'2011-11-08 14:54:56.163000','Munich','',3,'',47,'2010-02-26 22:57:13.310000');
INSERT INTO "users_profile" VALUES('','',48,'9a29dd1bf3d8aac62621c451e46a8d2a',2,0,'','',0,'2011-05-17 23:54:53.350000','Salt Lake City Utah','',3,'',48,'2010-02-26 23:38:01.027000');
INSERT INTO "users_profile" VALUES('','',49,'fafdcd838d9847a188e169b38ebb78b3',2,0,'','',0,'2010-02-27 00:20:47.790000','','',3,'',49,'2010-02-27 00:20:47.790000');
INSERT INTO "users_profile" VALUES('','',50,'cb1e8828b7d32e74d6809383d7639ec2',2,0,'','',0,'2010-03-01 23:25:38.710000','','',3,'',50,'2010-02-27 00:38:07.187000');
INSERT INTO "users_profile" VALUES('http://pbeltrao.blogspot.com','',51,'22fd820c589e1cae58caa9de512421d6',2,0,'','',0,'2010-04-21 02:16:04.980000','','',3,'',51,'2010-02-27 00:51:02.197000');
INSERT INTO "users_profile" VALUES('','',52,'58bdf64561f2558b9dbd9f5db93201a1',2,0,'','',0,'2011-06-27 13:20:27.173000','Johannesburg, South Africa','',3,'',52,'2010-02-27 01:00:28.793000');
INSERT INTO "users_profile" VALUES('http://www.michaelbarton.me.uk','<p>Doing post doc in genomics.</p>
',53,'66b43395f7793dbd8cc0d72ac942d401',2,0,'','',0,'2011-11-08 01:33:32.650000','Akron, Ohio, United States','',3,'',53,'2010-02-27 01:59:50.253000');
INSERT INTO "users_profile" VALUES('','',54,'e63472d230a169eac4b20ab349bab9a7',2,0,'','',0,'2011-11-07 16:47:03.237000','','',3,'',54,'2010-02-28 01:19:21.037000');
INSERT INTO "users_profile" VALUES('http://esysbio.bccs.uib.no/','<p>Just wanted to get the autobiographer badge ;)</p>
',55,'95f482a1800be96d221414dc23d607e8',2,0,'','',0,'2011-11-08 19:25:04.193000','Bergen','',3,'',55,'2010-02-28 03:00:55.783000');
INSERT INTO "users_profile" VALUES('http://bcbio.wordpress.com','',56,'1cb6356fb69fa99d9a66e723a201bc7e',2,0,'','',0,'2011-11-08 20:25:40.230000','Boston, MA','',3,'',56,'2010-03-01 04:41:45.450000');
INSERT INTO "users_profile" VALUES('','',57,'8914a3a1501821306d8b7417d3d82376',2,0,'','',0,'2010-03-01 19:31:45.800000','Buenos Aires, Argentina','',3,'',57,'2010-03-01 19:31:45.800000');
INSERT INTO "users_profile" VALUES('http://davebridges.github.com','',58,'eff9a70e718b56dee3553a33eb854b54',2,0,'','',0,'2011-11-08 20:16:50.120000','Ann Arbor, MI','',3,'',58,'2010-03-01 22:58:57.787000');
INSERT INTO "users_profile" VALUES('http://eridanus.net','<p>Senior NGS Computational Biologist at Oxford Gene Technology.</p><p></p>

<p>You can find out more about OGT''s products and services at <a href="http://www.ogt.co.uk" rel="nofollow">www.ogt.co.uk</a>
</p><p>
OGT offers a broad range of microarray and targeted sequencing services enabling high-quality, high-throughput genomic studies. The company also specialises in cytogenetic microarays and biomarker discovery and validation.
</p><p>
Blog: <a href="http://eridanus.net" rel="nofollow">eridanus.net</a></p><p>
Twitter: <a href="http://twitter.com/d_swan" rel="nofollow">twitter.com/d_swan</a>
</p><p>
I contribute to BioStar on my own time and therefore comments and opinions are my own and not that of my employer.</p>
',59,'6bc3aa8ad96c74f1fa56408aaa5aade2',2,0,'','',0,'2011-11-08 17:49:49.013000','Oxford, UK','',3,'',59,'2010-03-01 23:03:24.703000');
INSERT INTO "users_profile" VALUES('','',60,'0a0b2e0fb8d91f6193c201225111475b',2,0,'','',0,'2011-04-23 18:38:38.710000','Austria','',3,'',60,'2010-03-02 01:18:00.527000');
INSERT INTO "users_profile" VALUES('http://bioperl.org','<p>I''m very sleepy.</p>
',61,'45f4a41937c7c75c513ec455d9c85be0',2,0,'','',0,'2011-08-16 01:18:52.720000','University of Illinois Urbana-Champaign','',3,'',61,'2010-03-02 11:11:16.723000');
INSERT INTO "users_profile" VALUES('http://openwetware.org/wiki/User:Darek_Kedra','',62,'f8fc1be51e6abdb8b0c2456d6df6f048',2,0,'','',0,'2011-11-08 17:44:29.277000','Barcelona, Spain','',3,'',62,'2010-03-02 16:25:56.590000');
INSERT INTO "users_profile" VALUES('','',63,'bb238598571d678088579c7563df1e98',2,0,'','',0,'2011-03-18 16:15:05.320000','','',3,'',63,'2010-03-03 06:34:39.337000');
INSERT INTO "users_profile" VALUES('http://manuelcorpas.com','<p>I am a lead developer of the DECIPHER database, a “DatabasE of Chromosomal Imbalance and Phenotype in Humans using Ensembl Resources”. The DECIPHER database is hosted at the Wellcome Trust Sanger Institute. Its main aim is to collect genotype and phenotype patient information from hospitals and research institutions around the world to help diagnose and catalogue developmental diseases. I was awarded a PhD in Computer Science by the University of Manchester for an algorithm I developed to enhance the classification of protein families using sequence and structural signatures.</p>
',64,'9d62f6ddede1ab39b7ca09949c4003fc',2,0,'','',0,'2010-03-30 21:53:26.190000','Cambridge','',3,'',64,'2010-03-03 23:17:23.793000');
INSERT INTO "users_profile" VALUES('http://www.abhishek-tiwari.com','',65,'d465f52b3f36e823cfc5707da7c4bec7',2,0,'','',0,'2010-03-24 14:41:17.853000','','',3,'',65,'2010-03-04 03:18:31.757000');
INSERT INTO "users_profile" VALUES('http://nsaunders.wordpress.com','<p>I''m a statistical bioinformatician with CSIRO Mathematics, Information and Statistics (CMIS), based in Sydney, Australia.</p>
',66,'83f87d8f18e250d84cd68fcc5af296ad',2,0,'','',0,'2011-11-08 16:47:18.240000','Sydney, Australia','',3,'',66,'2010-03-04 13:55:19.360000');
INSERT INTO "users_profile" VALUES('http://friendfeed.com/byzia','',67,'20de3ac83169e18dbbb13115520a8f0a',2,0,'','',0,'2011-11-08 17:34:53.170000','Wroclaw, Poland','',3,'',67,'2010-03-04 15:00:00.847000');
INSERT INTO "users_profile" VALUES('http://jeroen.vangoey.be','<p>Just Another Genome Hacker, 
<br>
<br>
Bioinformatics software developer at <a href="http://www.applied-maths.com/" rel="nofollow">Applied Maths</a>
<br>
<br>
Also known as BioGeek on the web. (Twitter: <a href="http://twitter.com/BioGeek" rel="nofollow">@BioGeek</a>)</p>
',68,'6081840f76c59f8ae7ec4cdbde9e6b67',2,0,'','',0,'2011-11-08 02:04:57.533000','Antwerp, Belgium','',3,'',68,'2010-03-04 15:23:54.510000');
INSERT INTO "users_profile" VALUES('http://bioinformatics.ucdavis.edu','',69,'e96cc54bdc8be8a9ccac7e3a84f97716',2,0,'','',0,'2011-11-08 03:49:17.890000','Davis, CA','',3,'',69,'2010-03-04 15:49:54.847000');
INSERT INTO "users_profile" VALUES('http://friendfeed.com/walshtp','<p>Bioinformatics programmer and sysadmin. Increasingly grizzled old Perl hacker. Machine for turning tea into grumpiness.
</p>

<p>
Any views expressed here are mine alone and not necessarily representative of my employer.
</p>
',70,'875b9a284146ee3209a6859330cf55cd',2,0,'','',0,'2011-11-08 16:50:21.707000','','',3,'',70,'2010-03-04 18:26:30.240000');
INSERT INTO "users_profile" VALUES('http://egonw.github.com','<p>Post-doc working on Open Data, Open Source, and Open Standards in cheminformatics, chemometrics and bioinformatics. Projects I work on include the CDK, Bioclipse, the Blue Obelisk Data Repository, and much more…</p>
',71,'c262f57e4ffc3cf04195b45ae7ffc8ae',2,0,'','',0,'2011-11-08 19:30:22.413000','Uppsala','',3,'',71,'2010-03-04 21:17:25.543000');
INSERT INTO "users_profile" VALUES('http://mndoci.com','<p>Principal Product Manager for Amazon EC2.  Blogger, podcaster, scientist, technologist, musician and all round geek.  Find me at</p>

<p>

<a href="http://mndoci.com" rel="nofollow">bbgm</a> - the blog<br>
<a href="http://mndoci.github.com" rel="nofollow">mndoci.github.com</a> - the portal<br>
<a href="http://c2cbio.com" rel="nofollow">Coast to Coast Bio</a> - the podcast
</p>
',72,'6877dd5dd762e060565fac406bd27f07',2,0,'','',0,'2011-11-08 11:17:03.400000','','',3,'',72,'2010-03-04 21:17:55.313000');
INSERT INTO "users_profile" VALUES('http://jermdemo.blogspot.com','<p>Bioinformatics software developer working at a pediatric hospital.</p>
',73,'2c124c5550fe95c67d17a45d1d7a7c56',2,0,'','',0,'2011-11-08 20:01:40.577000','Philadelphia, PA','',3,'',73,'2010-03-04 21:31:15.363000');
INSERT INTO "users_profile" VALUES('','',74,'70f7b9a66ae472b5ed1b17f40d989eb7',2,0,'','',0,'2011-11-08 02:03:19.973000','','',3,'',74,'2010-03-04 21:35:18.017000');
INSERT INTO "users_profile" VALUES('http://onertipaday.blogspot.com/','',75,'cbd6bf00c3591e093928a98af1116ece',2,0,'','',0,'2011-11-08 16:22:19.363000','Trieste','',3,'',75,'2010-03-04 22:33:40.567000');
INSERT INTO "users_profile" VALUES('','',76,'06778e6e4d6ab6a67c9c9b53a6379cb0',2,0,'','',0,'2011-09-26 20:55:47.580000','','',3,'',76,'2010-03-04 22:42:20.013000');
INSERT INTO "users_profile" VALUES('http://eleanorhowe.com','',77,'8014b9d7bd420ac00764f9e71522d5ca',2,0,'','',0,'2010-05-14 02:55:19.540000','Boston','',3,'',77,'2010-03-04 22:53:55.957000');
INSERT INTO "users_profile" VALUES('http://bioinformatics.ucdavis.edu','',78,'c4fb8a744f8c1f98df65f6035481d954',2,0,'','',0,'2010-07-21 05:56:49.750000','Davis, CA','',3,'',78,'2010-03-05 00:05:38.647000');
INSERT INTO "users_profile" VALUES('','',79,'972a228cf87d1880992a4e49747934b3',2,0,'','',0,'2011-11-02 03:47:03.033000','Nijmegen, the Netherlands','',3,'',79,'2010-03-05 00:08:31.377000');
INSERT INTO "users_profile" VALUES('','',80,'0f28745fc2ff36f03b334997b2beb343',2,0,'','',0,'2011-11-08 20:07:57.847000','Vienna','',3,'',80,'2010-03-05 00:30:00.850000');
INSERT INTO "users_profile" VALUES('','',81,'0e14b6caad8e33ffd5e198fb79d6106f',2,0,'','',0,'2010-03-23 00:28:10.323000','','',3,'',81,'2010-03-05 01:06:35.720000');
INSERT INTO "users_profile" VALUES('','',82,'f6f43ac9a7c1011fe020bf3a8b3f3e79',2,0,'','',0,'2011-06-17 00:05:52.500000','MA','',3,'',82,'2010-03-05 01:56:05.877000');
INSERT INTO "users_profile" VALUES('','',83,'f726efb22afbb4736b604618bbaa50d1',2,0,'','',0,'2011-11-08 20:36:17.927000','','',3,'',83,'2010-03-05 02:11:58.703000');
INSERT INTO "users_profile" VALUES('','',84,'f8068b82cd9e489230bfd103213fb0ac',2,0,'','',0,'2010-03-12 21:52:48.487000','','',3,'',84,'2010-03-05 02:36:50.333000');
INSERT INTO "users_profile" VALUES('','',85,'f5b85830da966bcf4f0f1c0151bbd38d',2,0,'','',0,'2010-03-05 02:49:12.787000','','',3,'',85,'2010-03-05 02:49:12.787000');
INSERT INTO "users_profile" VALUES('http://sulab.org','<p>I am an Associate Professor at the Scripps Research Institute.  My lab is equally focused on biomedical discovery and bioinformatics tool development.</p>
',86,'bda0d772b394dd85d0742cc87ad6df5a',2,0,'','',0,'2011-11-08 20:18:14.087000','San Diego, CA','',3,'',86,'2010-03-05 05:27:54.383000');
INSERT INTO "users_profile" VALUES('http://twitter.com/kshameer','<p>Post-doc (Translational Bioinformatics) at Mayo Clinic, working on projects in computational genomics, network analysis and biomedical informatics. </p>

<p>PhD in Computational Biology (National Centre for Biological Sciences (NCBS - TIFR), Bangalore - India) </p>
',87,'a2e3b308bd3866d629083189b9a88929',2,0,'','',0,'2011-11-08 19:35:50.253000','Bangalore','',3,'',87,'2010-03-05 08:42:56.367000');
INSERT INTO "users_profile" VALUES('','',88,'494f41a8f4ae2a1654d067a05d0e11fe',2,0,'','',0,'2011-09-06 18:37:42.953000','','',3,'',88,'2010-03-05 21:59:48.507000');
INSERT INTO "users_profile" VALUES('http://stackoverflow.com/users/163080','<p>Bioinformatics research of brain tumor, mostly with MATLAB. </p>
',89,'6987baa9e233a80bea470fa56df7f1f1',2,0,'','',0,'2011-11-08 04:26:39.040000','Bethesda, MD','',3,'',89,'2010-03-06 07:42:33.603000');
INSERT INTO "users_profile" VALUES('','',90,'748035f1fd58d4cc501af44d8eda11d9',2,0,'','',0,'2011-10-14 08:50:24.157000','','',3,'',90,'2010-03-06 11:12:19.247000');
INSERT INTO "users_profile" VALUES('http://followthedata.wordpress.com','<p>twitter.com/mikaelhuss</p>
',91,'65590dcbf0a8db1225f8e957594cd9eb',2,0,'','',0,'2011-11-08 04:17:24.010000','Stockholm','',3,'',91,'2010-03-06 20:53:00.370000');
INSERT INTO "users_profile" VALUES('http://hostpathogen.blogspot.com/','<p>genomics &amp; computational biology student @ Upenn</p>
',92,'5b081e155365f6672e117f84556fb4a9',2,0,'','',0,'2011-11-02 01:42:19.240000','philadelphia','',3,'',92,'2010-03-06 21:35:50.630000');
INSERT INTO "users_profile" VALUES('','',93,'23394d913db5af989f0f3486b48e9081',2,0,'','',0,'2010-03-06 22:01:05.290000','','',3,'',93,'2010-03-06 22:01:05.290000');
INSERT INTO "users_profile" VALUES('http://rosettadesigngroup.com/blog/','',94,'4ccceae86eb24fc2c5b37320c7d32fe0',2,0,'','',0,'2010-06-13 15:58:45.870000','','',3,'',94,'2010-03-06 23:26:09.143000');
INSERT INTO "users_profile" VALUES('http://in.tegr.in','<p>Proven and highly competent engineer, biotechnology manager, and entrepreneur with over 8 years experience in providing intuitive project management, professional leadership, and deep technical understanding on a wide variety of biotechnology and web based projects.</p>
',95,'0b2b43d8271a1da6584bfc4803b5c7f4',2,0,'','',0,'2011-10-07 00:30:25.400000','DC Metro Area','',3,'',95,'2010-03-06 23:49:08.633000');
INSERT INTO "users_profile" VALUES('','',96,'d3e67a927b2ad89ca1385868756421f2',2,0,'','',0,'2010-04-15 00:10:50.287000','','',3,'',96,'2010-03-07 00:07:27.297000');
INSERT INTO "users_profile" VALUES('http://noble.gs.washington.edu/~mmh1/','',97,'9769659781d9111ddc2deaa67fcf4027',2,0,'','',0,'2010-04-28 03:33:36.557000','','',3,'',97,'2010-03-07 00:14:15.883000');
INSERT INTO "users_profile" VALUES('','',98,'1299ef07318c56f99580ad79edfeb8fe',2,0,'','',0,'2010-03-13 22:41:08.860000','','',3,'',98,'2010-03-07 00:20:20.613000');
INSERT INTO "users_profile" VALUES('http://sites.google.com/site/avilella','<p>I am at the Vertebrate Genomics team at EMBL-EBI in Cambridge, UK.</p>
',99,'2826f262603fe76fe05c1d461c1667ec',2,0,'','',0,'2011-11-08 16:10:20.983000','Cambridge, UK','',3,'',99,'2010-03-07 00:20:28.960000');
INSERT INTO "users_profile" VALUES('','',100,'c7d6b57e3851503cf47bafe0d32f82ea',2,0,'','',0,'2010-03-07 01:03:07.767000','','',3,'',100,'2010-03-07 01:03:07.767000');
INSERT INTO "users_profile" VALUES('','',101,'890350a827fe46c378b6f8ecc7d08288',2,0,'','',0,'2014-05-09 14:39:23.925000','','',3,'',101,'2014-05-09 14:39:23.887000');
CREATE TABLE "posts_tag" ("id" integer NOT NULL PRIMARY KEY, "name" text NOT NULL, "count" integer NOT NULL);
INSERT INTO "posts_tag" VALUES(1,'guidelines',7);
INSERT INTO "posts_tag" VALUES(2,'bed',14);
INSERT INTO "posts_tag" VALUES(3,'gff',7);
INSERT INTO "posts_tag" VALUES(4,'galaxy',18);
INSERT INTO "posts_tag" VALUES(5,'yeast',7);
INSERT INTO "posts_tag" VALUES(6,'motif',18);
INSERT INTO "posts_tag" VALUES(7,'microarray',7);
INSERT INTO "posts_tag" VALUES(8,'clustering',7);
INSERT INTO "posts_tag" VALUES(9,'test',7);
INSERT INTO "posts_tag" VALUES(10,'nucleotides',7);
INSERT INTO "posts_tag" VALUES(11,'solid',21);
INSERT INTO "posts_tag" VALUES(12,'deep',7);
INSERT INTO "posts_tag" VALUES(13,'sequence',42);
INSERT INTO "posts_tag" VALUES(14,'boy',14);
INSERT INTO "posts_tag" VALUES(15,'george',14);
INSERT INTO "posts_tag" VALUES(16,'geneid',7);
INSERT INTO "posts_tag" VALUES(17,'accession',7);
INSERT INTO "posts_tag" VALUES(18,'mapping',7);
INSERT INTO "posts_tag" VALUES(19,'conversion',14);
INSERT INTO "posts_tag" VALUES(20,'shrimp',7);
INSERT INTO "posts_tag" VALUES(21,'sequencing',7);
INSERT INTO "posts_tag" VALUES(22,'short',7);
INSERT INTO "posts_tag" VALUES(23,'aligner',7);
INSERT INTO "posts_tag" VALUES(24,'meme',7);
INSERT INTO "posts_tag" VALUES(25,'sge',7);
INSERT INTO "posts_tag" VALUES(26,'compilation',14);
INSERT INTO "posts_tag" VALUES(27,'rna',7);
INSERT INTO "posts_tag" VALUES(28,'general',28);
INSERT INTO "posts_tag" VALUES(29,'subjective',35);
INSERT INTO "posts_tag" VALUES(30,'os',7);
INSERT INTO "posts_tag" VALUES(31,'programming',7);
INSERT INTO "posts_tag" VALUES(32,'languages',7);
INSERT INTO "posts_tag" VALUES(33,'gene',7);
INSERT INTO "posts_tag" VALUES(34,'go',7);
INSERT INTO "posts_tag" VALUES(35,'string',7);
INSERT INTO "posts_tag" VALUES(36,'protein',14);
INSERT INTO "posts_tag" VALUES(37,'interacti',7);
INSERT INTO "posts_tag" VALUES(38,'ppi',7);
INSERT INTO "posts_tag" VALUES(39,'pin',7);
INSERT INTO "posts_tag" VALUES(40,'structure',7);
INSERT INTO "posts_tag" VALUES(41,'blast',7);
INSERT INTO "posts_tag" VALUES(42,'ucsc',7);
INSERT INTO "posts_tag" VALUES(43,'fasta',7);
INSERT INTO "posts_tag" VALUES(44,'hdf',7);
INSERT INTO "posts_tag" VALUES(45,'biohdf',7);
INSERT INTO "posts_tag" VALUES(46,'storage',7);
INSERT INTO "posts_tag" VALUES(47,'taverna',7);
INSERT INTO "posts_tag" VALUES(48,'plugin',7);
INSERT INTO "posts_tag" VALUES(49,'maven',7);
INSERT INTO "posts_tag" VALUES(50,'workflow',7);
INSERT INTO "posts_tag" VALUES(51,'format',7);
INSERT INTO "posts_tag" VALUES(52,'make',7);
INSERT INTO "posts_tag" VALUES(53,'pipeline',7);
INSERT INTO "posts_tag" VALUES(54,'organization',7);
INSERT INTO "posts_tag" VALUES(55,'agile',7);
INSERT INTO "posts_tag" VALUES(56,'good',7);
INSERT INTO "posts_tag" VALUES(57,'team',7);
INSERT INTO "posts_tag" VALUES(58,'biopython',7);
INSERT INTO "posts_tag" VALUES(59,'python',14);
INSERT INTO "posts_tag" VALUES(60,'pygr',7);
INSERT INTO "posts_tag" VALUES(61,'use',14);
INSERT INTO "posts_tag" VALUES(62,'interval',7);
INSERT INTO "posts_tag" VALUES(63,'query',7);
INSERT INTO "posts_tag" VALUES(64,'genomics',7);
CREATE TABLE "posts_post_tag_set" ("id" integer NOT NULL PRIMARY KEY, "post_id" integer NOT NULL, "tag_id" integer NOT NULL);
INSERT INTO "posts_post_tag_set" VALUES(1,1,1);
INSERT INTO "posts_post_tag_set" VALUES(2,2,3);
INSERT INTO "posts_post_tag_set" VALUES(3,2,2);
INSERT INTO "posts_post_tag_set" VALUES(4,2,4);
INSERT INTO "posts_post_tag_set" VALUES(5,4,5);
INSERT INTO "posts_post_tag_set" VALUES(6,4,6);
INSERT INTO "posts_post_tag_set" VALUES(7,5,8);
INSERT INTO "posts_post_tag_set" VALUES(8,5,7);
INSERT INTO "posts_post_tag_set" VALUES(9,6,9);
INSERT INTO "posts_post_tag_set" VALUES(10,10,10);
INSERT INTO "posts_post_tag_set" VALUES(11,13,11);
INSERT INTO "posts_post_tag_set" VALUES(12,13,13);
INSERT INTO "posts_post_tag_set" VALUES(13,13,12);
INSERT INTO "posts_post_tag_set" VALUES(14,20,15);
INSERT INTO "posts_post_tag_set" VALUES(15,20,14);
INSERT INTO "posts_post_tag_set" VALUES(16,21,15);
INSERT INTO "posts_post_tag_set" VALUES(17,21,14);
INSERT INTO "posts_post_tag_set" VALUES(18,22,19);
INSERT INTO "posts_post_tag_set" VALUES(19,22,18);
INSERT INTO "posts_post_tag_set" VALUES(20,22,17);
INSERT INTO "posts_post_tag_set" VALUES(21,22,16);
INSERT INTO "posts_post_tag_set" VALUES(22,24,20);
INSERT INTO "posts_post_tag_set" VALUES(23,24,21);
INSERT INTO "posts_post_tag_set" VALUES(24,24,22);
INSERT INTO "posts_post_tag_set" VALUES(25,24,23);
INSERT INTO "posts_post_tag_set" VALUES(26,28,24);
INSERT INTO "posts_post_tag_set" VALUES(27,28,25);
INSERT INTO "posts_post_tag_set" VALUES(28,28,26);
INSERT INTO "posts_post_tag_set" VALUES(29,28,6);
INSERT INTO "posts_post_tag_set" VALUES(30,31,11);
INSERT INTO "posts_post_tag_set" VALUES(31,31,27);
INSERT INTO "posts_post_tag_set" VALUES(32,31,4);
INSERT INTO "posts_post_tag_set" VALUES(33,33,30);
INSERT INTO "posts_post_tag_set" VALUES(34,33,28);
INSERT INTO "posts_post_tag_set" VALUES(35,33,29);
INSERT INTO "posts_post_tag_set" VALUES(36,34,32);
INSERT INTO "posts_post_tag_set" VALUES(37,34,31);
INSERT INTO "posts_post_tag_set" VALUES(38,34,29);
INSERT INTO "posts_post_tag_set" VALUES(39,41,33);
INSERT INTO "posts_post_tag_set" VALUES(40,41,29);
INSERT INTO "posts_post_tag_set" VALUES(41,41,34);
INSERT INTO "posts_post_tag_set" VALUES(42,46,39);
INSERT INTO "posts_post_tag_set" VALUES(43,46,38);
INSERT INTO "posts_post_tag_set" VALUES(44,46,37);
INSERT INTO "posts_post_tag_set" VALUES(45,46,29);
INSERT INTO "posts_post_tag_set" VALUES(46,46,35);
INSERT INTO "posts_post_tag_set" VALUES(47,46,36);
INSERT INTO "posts_post_tag_set" VALUES(48,48,13);
INSERT INTO "posts_post_tag_set" VALUES(49,48,36);
INSERT INTO "posts_post_tag_set" VALUES(50,48,40);
INSERT INTO "posts_post_tag_set" VALUES(51,51,13);
INSERT INTO "posts_post_tag_set" VALUES(52,51,41);
INSERT INTO "posts_post_tag_set" VALUES(53,53,11);
INSERT INTO "posts_post_tag_set" VALUES(54,56,13);
INSERT INTO "posts_post_tag_set" VALUES(55,56,42);
INSERT INTO "posts_post_tag_set" VALUES(56,56,43);
INSERT INTO "posts_post_tag_set" VALUES(57,58,28);
INSERT INTO "posts_post_tag_set" VALUES(58,58,29);
INSERT INTO "posts_post_tag_set" VALUES(59,69,46);
INSERT INTO "posts_post_tag_set" VALUES(60,69,44);
INSERT INTO "posts_post_tag_set" VALUES(61,69,45);
INSERT INTO "posts_post_tag_set" VALUES(62,76,47);
INSERT INTO "posts_post_tag_set" VALUES(63,76,26);
INSERT INTO "posts_post_tag_set" VALUES(64,76,50);
INSERT INTO "posts_post_tag_set" VALUES(65,76,48);
INSERT INTO "posts_post_tag_set" VALUES(66,76,49);
INSERT INTO "posts_post_tag_set" VALUES(67,77,19);
INSERT INTO "posts_post_tag_set" VALUES(68,77,2);
INSERT INTO "posts_post_tag_set" VALUES(69,77,51);
INSERT INTO "posts_post_tag_set" VALUES(70,79,54);
INSERT INTO "posts_post_tag_set" VALUES(71,79,28);
INSERT INTO "posts_post_tag_set" VALUES(72,79,53);
INSERT INTO "posts_post_tag_set" VALUES(73,79,52);
INSERT INTO "posts_post_tag_set" VALUES(74,88,55);
INSERT INTO "posts_post_tag_set" VALUES(75,88,57);
INSERT INTO "posts_post_tag_set" VALUES(76,88,56);
INSERT INTO "posts_post_tag_set" VALUES(77,88,28);
INSERT INTO "posts_post_tag_set" VALUES(78,90,13);
INSERT INTO "posts_post_tag_set" VALUES(79,90,59);
INSERT INTO "posts_post_tag_set" VALUES(80,90,58);
INSERT INTO "posts_post_tag_set" VALUES(81,92,60);
INSERT INTO "posts_post_tag_set" VALUES(82,92,61);
INSERT INTO "posts_post_tag_set" VALUES(83,92,13);
INSERT INTO "posts_post_tag_set" VALUES(84,92,59);
INSERT INTO "posts_post_tag_set" VALUES(85,99,61);
INSERT INTO "posts_post_tag_set" VALUES(86,99,62);
INSERT INTO "posts_post_tag_set" VALUES(87,99,63);
INSERT INTO "posts_post_tag_set" VALUES(88,99,64);
INSERT INTO "posts_post_tag_set" VALUES(89,101,4);
INSERT INTO "posts_post_tag_set" VALUES(90,101,6);
CREATE TABLE "posts_postview" ("id" integer NOT NULL PRIMARY KEY, "ip" char(39) NULL, "post_id" integer NOT NULL, "date" datetime NOT NULL);
INSERT INTO "posts_postview" VALUES(1,'127.0.0.1',34,'2014-04-29 17:28:09.157000');
INSERT INTO "posts_postview" VALUES(2,'127.0.0.1',92,'2014-04-29 17:28:09.779000');
INSERT INTO "posts_postview" VALUES(3,'127.0.0.1',34,'2014-04-30 17:10:36.535000');
INSERT INTO "posts_postview" VALUES(4,'127.0.0.1',101,'2014-05-09 14:39:20.173000');
INSERT INTO "posts_postview" VALUES(5,'127.0.0.1',99,'2014-05-09 14:40:31.856000');
INSERT INTO "posts_postview" VALUES(6,'127.0.0.1',69,'2014-07-11 14:15:13.590000');
INSERT INTO "posts_postview" VALUES(7,'127.0.0.1',99,'2014-07-11 14:16:28.239000');
CREATE TABLE "posts_vote" ("id" integer NOT NULL PRIMARY KEY, "author_id" integer NOT NULL, "post_id" integer NOT NULL, "type" integer NOT NULL, "date" datetime NOT NULL);
INSERT INTO "posts_vote" VALUES(1,2,7,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(2,2,7,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(3,2,4,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(4,2,8,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(5,10,11,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(6,10,11,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(7,10,11,3,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(8,10,12,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(9,2,10,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(10,2,13,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(11,2,15,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(12,2,13,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(13,2,16,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(14,10,1,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(15,2,18,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(16,2,18,3,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(17,5,10,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(18,2,24,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(19,2,25,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(20,14,25,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(21,14,26,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(22,2,27,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(23,2,28,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(24,20,29,3,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(25,2,30,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(26,2,34,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(27,2,33,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(28,2,31,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(29,2,39,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(30,23,8,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(31,23,7,0,'2014-04-29 15:02:17.740000');
INSERT INTO "posts_vote" VALUES(32,23,35,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(33,23,36,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(34,23,38,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(35,23,1,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(36,23,17,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(37,23,3,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(38,2,43,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(39,2,40,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(40,23,44,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(41,23,45,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(42,23,47,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(43,23,47,3,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(44,2,46,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(45,2,41,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(46,2,22,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(47,23,45,3,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(48,4,32,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(49,4,32,3,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(50,4,8,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(51,2,48,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(52,2,50,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(53,14,52,3,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(54,2,51,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(55,23,50,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(56,23,49,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(57,23,50,3,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(58,23,52,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(59,23,54,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(60,23,55,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(61,23,53,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(62,2,56,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(63,23,57,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(64,23,61,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(65,2,59,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(66,2,56,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(67,2,64,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(68,2,38,1,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(69,2,63,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(70,2,65,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(71,2,66,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(72,33,46,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(73,33,60,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(74,33,4,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(75,31,62,0,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(76,31,2,2,'2014-04-29 15:02:17.741000');
INSERT INTO "posts_vote" VALUES(77,33,34,0,'2014-04-29 15:02:17.742000');
INSERT INTO "posts_vote" VALUES(78,33,56,0,'2014-04-29 15:02:17.742000');
INSERT INTO "posts_vote" VALUES(79,23,62,0,'2014-04-29 15:02:17.742000');
INSERT INTO "posts_vote" VALUES(80,23,64,0,'2014-04-29 15:02:17.742000');
INSERT INTO "posts_vote" VALUES(81,23,59,3,'2014-04-29 15:02:17.742000');
INSERT INTO "posts_vote" VALUES(82,23,59,0,'2014-04-29 15:02:17.742000');
INSERT INTO "posts_vote" VALUES(83,23,67,0,'2014-04-29 15:02:17.742000');
INSERT INTO "posts_vote" VALUES(84,23,63,0,'2014-04-29 15:02:17.742000');
INSERT INTO "posts_vote" VALUES(85,23,63,3,'2014-04-29 15:02:17.742000');
INSERT INTO "posts_vote" VALUES(86,23,27,0,'2014-04-29 15:02:17.742000');
INSERT INTO "posts_vote" VALUES(87,23,65,0,'2014-04-29 15:02:17.742000');
INSERT INTO "posts_vote" VALUES(88,23,68,0,'2014-04-29 15:02:17.742000');
INSERT INTO "posts_vote" VALUES(89,38,58,0,'2014-04-29 15:02:17.742000');
INSERT INTO "posts_vote" VALUES(90,23,68,3,'2014-04-29 15:02:17.742000');
INSERT INTO "posts_vote" VALUES(91,30,68,0,'2014-04-29 15:02:17.742000');
INSERT INTO "posts_vote" VALUES(92,23,69,0,'2014-04-29 15:02:17.742000');
INSERT INTO "posts_vote" VALUES(93,23,70,0,'2014-04-29 15:02:17.742000');
INSERT INTO "posts_vote" VALUES(94,2,67,0,'2014-04-29 15:02:17.742000');
INSERT INTO "posts_vote" VALUES(95,2,70,0,'2014-04-29 15:02:17.742000');
INSERT INTO "posts_vote" VALUES(96,2,71,0,'2014-04-29 15:02:17.742000');
INSERT INTO "posts_vote" VALUES(97,23,71,0,'2014-04-29 15:02:17.742000');
INSERT INTO "posts_vote" VALUES(98,23,73,0,'2014-04-29 15:02:17.742000');
INSERT INTO "posts_vote" VALUES(99,2,69,0,'2014-04-29 15:02:17.742000');
INSERT INTO "posts_vote" VALUES(100,30,73,0,'2014-04-29 15:02:17.742000');
CREATE TABLE "posts_subscription" ("id" integer NOT NULL PRIMARY KEY, "user_id" integer NOT NULL, "post_id" integer NOT NULL, "type" integer NOT NULL, "date" datetime NOT NULL);
INSERT INTO "posts_subscription" VALUES(1,2,1,1,'2014-04-29 15:02:13.939000');
INSERT INTO "posts_subscription" VALUES(2,2,2,1,'2014-04-29 15:02:13.973000');
INSERT INTO "posts_subscription" VALUES(3,3,4,1,'2014-04-29 15:02:14.052000');
INSERT INTO "posts_subscription" VALUES(4,2,5,1,'2014-04-29 15:02:14.084000');
INSERT INTO "posts_subscription" VALUES(5,5,6,1,'2014-04-29 15:02:14.118000');
INSERT INTO "posts_subscription" VALUES(6,5,4,0,'2014-04-29 15:02:14.155000');
INSERT INTO "posts_subscription" VALUES(7,6,4,0,'2014-04-29 15:02:14.184000');
INSERT INTO "posts_subscription" VALUES(8,7,4,0,'2014-04-29 15:02:14.211000');
INSERT INTO "posts_subscription" VALUES(9,10,10,1,'2014-04-29 15:02:14.234000');
INSERT INTO "posts_subscription" VALUES(10,2,10,0,'2014-04-29 15:02:14.280000');
INSERT INTO "posts_subscription" VALUES(11,12,13,1,'2014-04-29 15:02:14.330000');
INSERT INTO "posts_subscription" VALUES(12,2,13,0,'2014-04-29 15:02:14.372000');
INSERT INTO "posts_subscription" VALUES(13,5,13,0,'2014-04-29 15:02:14.400000');
INSERT INTO "posts_subscription" VALUES(14,13,13,0,'2014-04-29 15:02:14.426000');
INSERT INTO "posts_subscription" VALUES(15,10,1,0,'2014-04-29 15:02:14.456000');
INSERT INTO "posts_subscription" VALUES(16,14,1,0,'2014-04-29 15:02:14.485000');
INSERT INTO "posts_subscription" VALUES(17,5,5,0,'2014-04-29 15:02:14.523000');
INSERT INTO "posts_subscription" VALUES(18,15,20,1,'2014-04-29 15:02:14.549000');
INSERT INTO "posts_subscription" VALUES(19,15,21,1,'2014-04-29 15:02:14.582000');
INSERT INTO "posts_subscription" VALUES(20,16,22,1,'2014-04-29 15:02:14.632000');
INSERT INTO "posts_subscription" VALUES(21,2,22,0,'2014-04-29 15:02:14.676000');
INSERT INTO "posts_subscription" VALUES(22,14,24,1,'2014-04-29 15:02:14.703000');
INSERT INTO "posts_subscription" VALUES(23,18,24,0,'2014-04-29 15:02:14.746000');
INSERT INTO "posts_subscription" VALUES(24,2,24,0,'2014-04-29 15:02:14.777000');
INSERT INTO "posts_subscription" VALUES(25,19,22,0,'2014-04-29 15:02:14.808000');
INSERT INTO "posts_subscription" VALUES(26,20,28,1,'2014-04-29 15:02:14.844000');
INSERT INTO "posts_subscription" VALUES(27,2,28,0,'2014-04-29 15:02:14.891000');
INSERT INTO "posts_subscription" VALUES(28,20,2,0,'2014-04-29 15:02:15.006000');
INSERT INTO "posts_subscription" VALUES(29,4,31,1,'2014-04-29 15:02:15.029000');
INSERT INTO "posts_subscription" VALUES(30,2,31,0,'2014-04-29 15:02:15.071000');
INSERT INTO "posts_subscription" VALUES(31,23,33,1,'2014-04-29 15:02:15.096000');
INSERT INTO "posts_subscription" VALUES(32,23,34,1,'2014-04-29 15:02:15.137000');
INSERT INTO "posts_subscription" VALUES(33,2,33,0,'2014-04-29 15:02:15.179000');
INSERT INTO "posts_subscription" VALUES(34,2,34,0,'2014-04-29 15:02:15.208000');
INSERT INTO "posts_subscription" VALUES(35,25,34,0,'2014-04-29 15:02:15.268000');
INSERT INTO "posts_subscription" VALUES(36,23,22,0,'2014-04-29 15:02:15.299000');
INSERT INTO "posts_subscription" VALUES(37,23,4,0,'2014-04-29 15:02:15.328000');
INSERT INTO "posts_subscription" VALUES(38,23,41,1,'2014-04-29 15:02:15.356000');
INSERT INTO "posts_subscription" VALUES(39,23,1,0,'2014-04-29 15:02:15.393000');
INSERT INTO "posts_subscription" VALUES(40,2,41,0,'2014-04-29 15:02:15.457000');
INSERT INTO "posts_subscription" VALUES(41,4,41,0,'2014-04-29 15:02:15.483000');
INSERT INTO "posts_subscription" VALUES(42,23,46,1,'2014-04-29 15:02:15.508000');
INSERT INTO "posts_subscription" VALUES(43,2,46,0,'2014-04-29 15:02:15.563000');
INSERT INTO "posts_subscription" VALUES(44,23,48,1,'2014-04-29 15:02:15.588000');
INSERT INTO "posts_subscription" VALUES(45,2,48,0,'2014-04-29 15:02:15.627000');
INSERT INTO "posts_subscription" VALUES(46,26,48,0,'2014-04-29 15:02:15.659000');
INSERT INTO "posts_subscription" VALUES(47,14,51,1,'2014-04-29 15:02:15.701000');
INSERT INTO "posts_subscription" VALUES(48,2,51,0,'2014-04-29 15:02:15.739000');
INSERT INTO "posts_subscription" VALUES(49,27,53,1,'2014-04-29 15:02:15.772000');
INSERT INTO "posts_subscription" VALUES(50,27,31,0,'2014-04-29 15:02:15.816000');
INSERT INTO "posts_subscription" VALUES(51,28,53,0,'2014-04-29 15:02:15.847000');
INSERT INTO "posts_subscription" VALUES(52,23,56,1,'2014-04-29 15:02:15.871000');
INSERT INTO "posts_subscription" VALUES(53,2,56,0,'2014-04-29 15:02:15.911000');
INSERT INTO "posts_subscription" VALUES(54,23,58,1,'2014-04-29 15:02:15.939000');
INSERT INTO "posts_subscription" VALUES(55,30,56,0,'2014-04-29 15:02:15.977000');
INSERT INTO "posts_subscription" VALUES(56,31,4,0,'2014-04-29 15:02:16.005000');
INSERT INTO "posts_subscription" VALUES(57,31,58,0,'2014-04-29 15:02:16.034000');
INSERT INTO "posts_subscription" VALUES(58,2,58,0,'2014-04-29 15:02:16.063000');
INSERT INTO "posts_subscription" VALUES(59,35,46,0,'2014-04-29 15:02:16.121000');
INSERT INTO "posts_subscription" VALUES(60,35,34,0,'2014-04-29 15:02:16.150000');
INSERT INTO "posts_subscription" VALUES(61,33,33,0,'2014-04-29 15:02:16.178000');
INSERT INTO "posts_subscription" VALUES(62,33,1,0,'2014-04-29 15:02:16.227000');
INSERT INTO "posts_subscription" VALUES(63,38,33,0,'2014-04-29 15:02:16.267000');
INSERT INTO "posts_subscription" VALUES(64,38,58,0,'2014-04-29 15:02:16.302000');
INSERT INTO "posts_subscription" VALUES(65,30,69,1,'2014-04-29 15:02:16.328000');
INSERT INTO "posts_subscription" VALUES(66,40,58,0,'2014-04-29 15:02:16.373000');
INSERT INTO "posts_subscription" VALUES(67,42,69,0,'2014-04-29 15:02:16.407000');
INSERT INTO "posts_subscription" VALUES(68,23,69,0,'2014-04-29 15:02:16.438000');
INSERT INTO "posts_subscription" VALUES(69,2,69,0,'2014-04-29 15:02:16.475000');
INSERT INTO "posts_subscription" VALUES(70,26,1,0,'2014-04-29 15:02:16.558000');
INSERT INTO "posts_subscription" VALUES(71,30,76,1,'2014-04-29 15:02:16.589000');
INSERT INTO "posts_subscription" VALUES(72,44,77,1,'2014-04-29 15:02:16.628000');
INSERT INTO "posts_subscription" VALUES(73,10,77,0,'2014-04-29 15:02:16.672000');
INSERT INTO "posts_subscription" VALUES(74,23,79,1,'2014-04-29 15:02:16.694000');
INSERT INTO "posts_subscription" VALUES(75,2,79,0,'2014-04-29 15:02:16.736000');
INSERT INTO "posts_subscription" VALUES(76,47,79,0,'2014-04-29 15:02:16.815000');
INSERT INTO "posts_subscription" VALUES(77,7,48,0,'2014-04-29 15:02:16.842000');
INSERT INTO "posts_subscription" VALUES(78,53,79,0,'2014-04-29 15:02:16.879000');
INSERT INTO "posts_subscription" VALUES(79,24,58,0,'2014-04-29 15:02:16.995000');
INSERT INTO "posts_subscription" VALUES(80,7,58,0,'2014-04-29 15:02:17.107000');
INSERT INTO "posts_subscription" VALUES(81,55,69,0,'2014-04-29 15:02:17.193000');
INSERT INTO "posts_subscription" VALUES(82,23,88,1,'2014-04-29 15:02:17.228000');
INSERT INTO "posts_subscription" VALUES(83,2,88,0,'2014-04-29 15:02:17.296000');
INSERT INTO "posts_subscription" VALUES(84,10,90,1,'2014-04-29 15:02:17.335000');
INSERT INTO "posts_subscription" VALUES(85,39,88,0,'2014-04-29 15:02:17.403000');
INSERT INTO "posts_subscription" VALUES(86,10,92,1,'2014-04-29 15:02:17.432000');
INSERT INTO "posts_subscription" VALUES(87,45,5,0,'2014-04-29 15:02:17.470000');
INSERT INTO "posts_subscription" VALUES(88,23,90,0,'2014-04-29 15:02:17.498000');
INSERT INTO "posts_subscription" VALUES(89,24,79,0,'2014-04-29 15:02:17.536000');
INSERT INTO "posts_subscription" VALUES(90,24,90,0,'2014-04-29 15:02:17.573000');
INSERT INTO "posts_subscription" VALUES(91,59,22,0,'2014-04-29 15:02:17.602000');
INSERT INTO "posts_subscription" VALUES(92,55,22,0,'2014-04-29 15:02:17.638000');
INSERT INTO "posts_subscription" VALUES(93,10,99,1,'2014-04-29 15:02:17.666000');
INSERT INTO "posts_subscription" VALUES(94,2,99,0,'2014-04-29 15:02:17.719000');
INSERT INTO "posts_subscription" VALUES(95,1,101,1,'2014-05-09 14:39:20.059000');
INSERT INTO "posts_subscription" VALUES(96,2,101,1,'2014-05-09 14:39:20.066000');
CREATE TABLE "posts_post" ("site_id" integer, "rank" real NOT NULL, "creation_date" datetime NOT NULL, "reply_count" integer NOT NULL, "tag_val" varchar(100) NOT NULL, "id" integer PRIMARY KEY, "view_count" integer NOT NULL, "thread_score" integer NOT NULL, "title" varchar(200) NOT NULL, "has_accepted" bool NOT NULL, "vote_count" integer NOT NULL, "content" text NOT NULL, "parent_id" integer, "comment_count" integer NOT NULL, "html" text NOT NULL, "type" integer NOT NULL, "status" integer NOT NULL, "book_count" integer NOT NULL, "root_id" integer, "lastedit_user_id" integer NOT NULL, "sticky" bool NOT NULL, "lastedit_date" datetime NOT NULL, "changed" bool NOT NULL, "subs_count" integer NOT NULL, "author_id" integer NOT NULL);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-09-30 20:12:07.053000',5,'guidelines',1,70,8,'Site Use Guidelines',0,2,'<p>Here are a few guidelines:</p>

<ol>
<li>The site''s goal is to answer bioinformatics and systems biology related questions</li>
<li>Answer questions to gain <em>reputation</em>. </li>
<li>Don''t forget to vote for answers that you like! Registered users may vote on answers.</li>
<li>If you are the one asking the original question you may also select the best answer</li>
<li>Subscribe to the RSS feeds for all questions or a single question to keep up to date with the developments</li>
</ol>
',1,0,'<p>Here are a few guidelines:</p>

<ol>
<li>The site''s goal is to answer bioinformatics and systems biology related questions</li>
<li>Answer questions to gain <em>reputation</em>. </li>
<li>Don''t forget to vote for answers that you like! Registered users may vote on answers.</li>
<li>If you are the one asking the original question you may also select the best answer</li>
<li>Subscribe to the RSS feeds for all questions or a single question to keep up to date with the developments</li>
</ol>
',0,1,0,1,2,0,'2010-02-26 20:10:59.777000',1,42,2);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-09-30 20:55:00.133000',2,'bed,gff,galaxy',2,187,3,'How Do I Convert From Bed Format To Gff Format?',0,0,'<p>I have a file in GFF format and I need to convert it to BED format. What do I do?</p>
',2,0,'<p>I have a file in GFF format and I need to convert it to BED format. What do I do?</p>
',0,1,1,2,2,0,'2009-09-30 21:09:43.357000',1,14,2);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-09-30 20:56:18.353000',0,'',3,0,0,'A: How Do I Convert From Bed Format To Gff Format?',0,1,'<p>Both formats are tab delimited text files used to represent DNA features in genomes. The order of columns between the two are different, there are also columns that correspond to attributes missing from one or the other format. Nonetheless <strong>the most important</strong> difference between the two is the coordinate systems that they assume. </p>

<p>The <a href=''http://genome.ucsc.edu/FAQ/FAQformat#format1''>BED</a> format developed at <code>UCSC</code> uses a zero based indexing and an open end interval whereas the <a href=''http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml''>GFF</a> format developed at <code>Sanger</code> assumes a 1 based coordinate system that includes both start and end coordinates. Therefore</p>

<p>The <code>[0,100]</code> interval in <code>BED</code> format corresponds to <code>[1,100]</code> in <code>GFF</code> format and both are <code>100</code> base long. That the first element in BED format will be have the index of <code>0</code> where the last <code>100th</code> element will have the index of <code>99</code>! Whereas in <code>GFF</code> the first element will have the index of <code>1</code> and the last element will have the index of <code>100</code>.</p>

<p>To convert between the two you may use <a href=''http://main.g2.bx.psu.edu/''>Galaxy</a> and select the section called <code>Select Formats</code> that will list various transformation options.</p>
',2,0,'<p>Both formats are tab delimited text files used to represent DNA features in genomes. The order of columns between the two are different, there are also columns that correspond to attributes missing from one or the other format. Nonetheless <strong>the most important</strong> difference between the two is the coordinate systems that they assume. </p>

<p>The <a rel="nofollow" href="http://genome.ucsc.edu/FAQ/FAQformat#format1">BED</a> format developed at <code>UCSC</code> uses a zero based indexing and an open end interval whereas the <a rel="nofollow" href="http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml">GFF</a> format developed at <code>Sanger</code> assumes a 1 based coordinate system that includes both start and end coordinates. Therefore</p>

<p>The <code>[0,100]</code> interval in <code>BED</code> format corresponds to <code>[1,100]</code> in <code>GFF</code> format and both are <code>100</code> base long. That the first element in BED format will be have the index of <code>0</code> where the last <code>100th</code> element will have the index of <code>99</code>! Whereas in <code>GFF</code> the first element will have the index of <code>1</code> and the last element will have the index of <code>100</code>.</p>

<p>To convert between the two you may use <a rel="nofollow" href="http://main.g2.bx.psu.edu/">Galaxy</a> and select the section called <code>Select Formats</code> that will list various transformation options.</p>
',1,1,0,2,2,0,'2009-10-01 01:46:25.510000',1,0,2);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-09-30 22:09:06.677000',5,'yeast,motif',4,150,11,'Finding Common Motifs In Sequences',0,3,'<p>I have a few hundred yeast sequences (20-80bp long) and I want to find common motifs (conserved bases at certain indices) in them. I am using a Mac</p>
',4,0,'<p>I have a few hundred yeast sequences (20-80bp long) and I want to find common motifs (conserved bases at certain indices) in them. I am using a Mac</p>
',0,1,0,4,3,0,'2009-09-30 22:09:06.677000',1,42,3);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-09-30 22:44:22.647000',3,'microarray,clustering',5,154,5,'Recommend Easy To Use Microarray Clustering Software',0,1,'<p>Feel free to post your favorite clustering tool.</p>
',5,0,'<p>Feel free to post your favorite clustering tool.</p>
',0,1,0,5,2,0,'2009-10-05 22:10:10.997000',1,21,2);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-10-01 00:49:39.563000',0,'test',6,2,0,'Test By Zhenhai',0,0,'<p>Hi, I just created my user id a few minutes ago. </p>

<p>Post this question to see how it works.</p>
',6,0,'<p>Hi, I just created my user id a few minutes ago. </p>

<p>Post this question to see how it works.</p>
',0,3,0,6,2,0,'2009-10-01 00:57:55.260000',1,7,5);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-10-01 00:55:02.457000',0,'',7,0,0,'A: Finding Common Motifs In Sequences',0,3,'<p>try this out?</p>

<p><a href=''http://fraenkel.mit.edu/webmotifs/form.html''>http://fraenkel.mit.edu/webmotifs/form.html</a></p>
',4,0,'<p>try this out?</p>

<p><a rel="nofollow" href="http://fraenkel.mit.edu/webmotifs/form.html">http://fraenkel.mit.edu/webmotifs/form.html</a></p>
',1,1,0,4,5,0,'2009-10-01 00:55:02.457000',1,0,5);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-10-01 01:32:29.097000',0,'',8,0,0,'A: Finding Common Motifs In Sequences',0,3,'<p>You can also use MEME:  <a href=''http://meme.sdsc.edu/''>http://meme.sdsc.edu/</a>.</p>
',4,0,'<p>You can also use MEME:  <a rel="nofollow" href="http://meme.sdsc.edu/">http://meme.sdsc.edu/</a>.</p>
',1,1,0,4,6,0,'2009-10-01 01:32:29.097000',1,0,6);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-10-01 01:35:28.020000',0,'',9,0,0,'A: Finding Common Motifs In Sequences',0,0,'<pre>
ACGGGCCCGACGATGCGTCGTA

ACGTACGTCGAACCGTCGTCGT

ACGTGCGTCGAAACGTCAGTCG

ACGGGTTCGATCGTCGTCGTCG
</pre>

<p>may be in Python I will break down the first sequence of required motif length into a sliding window and will search for those list of motifs in the rest of sequences using regular expression in python using <code>re.search()</code> method.</p>
',4,0,'<pre>ACGGGCCCGACGATGCGTCGTA

ACGTACGTCGAACCGTCGTCGT

ACGTGCGTCGAAACGTCAGTCG

ACGGGTTCGATCGTCGTCGTCG
</pre>

<p>may be in Python I will break down the first sequence of required motif length into a sliding window and will search for those list of motifs in the rest of sequences using regular expression in python using <code>re.search()</code> method.</p>
',1,1,0,4,2,0,'2009-10-01 07:09:51.160000',1,0,7);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-10-05 21:51:37.043000',1,'nucleotides',10,18,4,'How To Generate Multi-Nucleotide Occupancy Counts For Each Coordinate Of My Reads?',0,2,'<p>I need to generate nucleotide occupancy counts for each position of a given sequence then summed over each of the input sequences. An example desired output (for di-nucleotide AT):</p>

<p><img src=''http://github.com/ialbert/biostar-codesample/raw/master/python/images/dinuc.png'' alt=''dinucleotide occupancy'' /></p>
',10,0,'<p>I need to generate nucleotide occupancy counts for each position of a given sequence then summed over each of the input sequences. An example desired output (for di-nucleotide AT):</p>

<p><img src="http://github.com/ialbert/biostar-codesample/raw/master/python/images/dinuc.png" alt="dinucleotide occupancy"></p>
',0,1,0,10,2,0,'2009-10-13 20:53:15.190000',1,14,10);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-10-05 22:03:01.060000',0,'',11,0,0,'A: How To Generate Multi-Nucleotide Occupancy Counts For Each Coordinate Of My Read',0,2,'<p>The code snippet below will populate the <code>store</code> dictionary keyed by the nucleotide patterns and values as lists that contain the occupancy for each index. (Updated answer now includes arbitrary lenght nucleotide counts)::</p>

<pre><code>from itertools import count

def pattern_update(sequence, width=2, store={}):
    """
    Accumulates nucleotide patterns of a certain width with 
    position counts at each index.
    """

    # open intervals need a padding at end for proper slicing
    size  = len(sequence) + 1

    def zeroes():
        "Generates an empty array that holds the positions"
        return [ 0 ] * (size - width)

    # these are the end indices
    ends = range(width, size)

    for lo, hi in zip(count(), ends):
        # upon encoutering a missing key initialize 
        # that value for that key to the return value of the empty() function
        key = sequence[lo:hi]
        store.setdefault(key, zeroes())[lo] += 1

    return store
</code></pre>

<p>The code at <a href=''http://github.com/ialbert/biostar-codesample/blob/master/python/multipatt.py''>multipatt.py</a> demonstrates its use in a full program. Set the <code>size</code> to the maximal possible sequence size. A typical use case::</p>

<pre><code>store = {}
seq1 = ''ATGCT''
pattern_update(seq1, width=2, store=store)    

seq2 = ''ATCGC''
pattern_update(seq2, width=2, store=store)    

print store
</code></pre>

<p>will print::</p>

<pre><code>{''CG'': [0, 0, 1, 0], ''GC'': [0, 0, 1, 1], ''AT'': [2, 0, 0, 0], 
''TG'': [0, 1, 0, 0], ''TC'': [0, 1, 0, 0], ''CT'': [0, 0, 0, 1]}
</code></pre>
',10,0,'<p>The code snippet below will populate the <code>store</code> dictionary keyed by the nucleotide patterns and values as lists that contain the occupancy for each index. (Updated answer now includes arbitrary lenght nucleotide counts)::</p>

<pre><code>from itertools import count

def pattern_update(sequence, width=2, store={}):
    """
    Accumulates nucleotide patterns of a certain width with 
    position counts at each index.
    """

    # open intervals need a padding at end for proper slicing
    size  = len(sequence) + 1

    def zeroes():
        "Generates an empty array that holds the positions"
        return [ 0 ] * (size - width)

    # these are the end indices
    ends = range(width, size)

    for lo, hi in zip(count(), ends):
        # upon encoutering a missing key initialize 
        # that value for that key to the return value of the empty() function
        key = sequence[lo:hi]
        store.setdefault(key, zeroes())[lo] += 1

    return store
</code></pre>

<p>The code at <a rel="nofollow" href="http://github.com/ialbert/biostar-codesample/blob/master/python/multipatt.py">multipatt.py</a> demonstrates its use in a full program. Set the <code>size</code> to the maximal possible sequence size. A typical use case::</p>

<pre><code>store = {}
seq1 = ''ATGCT''
pattern_update(seq1, width=2, store=store)    

seq2 = ''ATCGC''
pattern_update(seq2, width=2, store=store)    

print store
</code></pre>

<p>will print::</p>

<pre><code>{''CG'': [0, 0, 1, 0], ''GC'': [0, 0, 1, 1], ''AT'': [2, 0, 0, 0], 
''TG'': [0, 1, 0, 0], ''TC'': [0, 1, 0, 0], ''CT'': [0, 0, 0, 1]}
</code></pre>
',1,1,0,10,2,0,'2009-10-13 21:00:00.047000',1,0,2);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-10-05 22:09:57.673000',0,'',12,0,0,'A: Recommend Easy To Use Microarray Clustering Software',0,2,'<p>One of my favorites is the <a href=''http://www.tm4.org/mev.html''>MEV</a> micro-array data analysis tool.
It is simple to use and it has a very large number of features. </p>

<p>Works well for any type of data. You can also load into it data from a file that is in a simple text format:</p>

<pre>
GENE1, value1, value2
GENE2, value1, value2
</pre>

<p>Feel free to post your favorite clustering tool.</p>
',5,0,'<p>One of my favorites is the <a rel="nofollow" href="http://www.tm4.org/mev.html">MEV</a> micro-array data analysis tool.
It is simple to use and it has a very large number of features. </p>

<p>Works well for any type of data. You can also load into it data from a file that is in a simple text format:</p>

<pre>GENE1, value1, value2
GENE2, value1, value2
</pre>

<p>Feel free to post your favorite clustering tool.</p>
',1,1,0,5,2,0,'2009-10-05 22:09:57.673000',1,0,2);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-10-07 00:58:10.227000',3,'solid,deep,sequence',13,72,4,'Chip Dna Deep Sequence',0,2,'<p>Hi, everyone,
I am posting this question for my friend.
He is analyzing his CHIP DNA solid deep sequence data, and find out that near 80% reads can not be mapped to the human genome. We are wondering if this high percentage unmapped reads is normal in CHIP DNA deep sequence or there may be something wrong with his result.</p>
',13,0,'<p>Hi, everyone,
I am posting this question for my friend.
He is analyzing his CHIP DNA solid deep sequence data, and find out that near 80% reads can not be mapped to the human genome. We are wondering if this high percentage unmapped reads is normal in CHIP DNA deep sequence or there may be something wrong with his result.</p>
',0,1,0,13,12,0,'2009-10-07 00:58:10.227000',1,28,12);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-10-07 02:10:32.730000',0,'',14,0,0,'A: Chip Dna Deep Sequence',0,0,'<p>I recall that our first samples that we ran on the Solid sequencer have had bad performance. Not quite an 80% loss but around 40%-60% reads were unmappable (yeast). Some other lab members will hopefully chime in with more details. </p>
',13,0,'<p>I recall that our first samples that we ran on the Solid sequencer have had bad performance. Not quite an 80% loss but around 40%-60% reads were unmappable (yeast). Some other lab members will hopefully chime in with more details. </p>
',1,1,0,13,2,0,'2009-10-07 02:10:32.730000',1,0,2);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-10-07 02:41:24.877000',0,'',15,0,0,'A: Chip Dna Deep Sequence',0,1,'<p>Hi there,</p>

<p>We have done numbers of SOLiD sequencing run on yeast samples. Normally there are only 30-40 percent of total tags can be uniquely mapped back to yeast genome. </p>

<p>What I would recommend is do it on solexa. You get much higher quality tags.</p>

<p>cheers,</p>
',13,0,'<p>Hi there,</p>

<p>We have done numbers of SOLiD sequencing run on yeast samples. Normally there are only 30-40 percent of total tags can be uniquely mapped back to yeast genome. </p>

<p>What I would recommend is do it on solexa. You get much higher quality tags.</p>

<p>cheers,</p>
',1,1,0,13,5,0,'2009-10-07 02:41:24.877000',1,0,5);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-10-07 04:04:41.747000',0,'',16,0,0,'A: Chip Dna Deep Sequence',0,1,'<p>Your 20% mapping yield looks like low for normal ChIP experiment, even for human. Several factors can reduce this mapping yield. I am wondering which kind of ChIP was used in your case. That is, which kind of proteins was ChIPed?</p>
',13,0,'<p>Your 20% mapping yield looks like low for normal ChIP experiment, even for human. Several factors can reduce this mapping yield. I am wondering which kind of ChIP was used in your case. That is, which kind of proteins was ChIPed?</p>
',1,1,0,13,13,0,'2009-10-07 04:04:41.747000',1,0,13);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-10-07 22:41:44.823000',0,'',17,0,0,'A: Site Use Guidelines',0,1,'<p>If you are shy about asking the question on your own behalf submit it to to the <strong>Question Bot</strong> and it will be posted anonymously. Send email to the <code>Question Bot</code> link at the bottom.</p>
',1,0,'<p>If you are shy about asking the question on your own behalf submit it to to the <strong>Question Bot</strong> and it will be posted anonymously. Send email to the <code>Question Bot</code> link at the bottom.</p>
',1,1,0,1,2,0,'2009-10-07 23:16:18.690000',1,0,10);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-10-09 09:28:20.413000',0,'',18,0,0,'A: Site Use Guidelines',0,1,'<p>Hi,</p>

<p>I don''t think a new user can vote on a question or an answer.
The site says I need 15 reputation...</p>
',1,0,'<p>Hi,</p>

<p>I don''t think a new user can vote on a question or an answer.
The site says I need 15 reputation...</p>
',1,1,0,1,14,0,'2009-10-09 09:28:20.413000',1,0,14);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-10-16 18:25:38.237000',0,'',19,0,0,'A: Recommend Easy To Use Microarray Clustering Software',0,1,'<p>I would recommend a combination of cluster and treeview.</p>

<p>pretty powerful!</p>
',5,0,'<p>I would recommend a combination of cluster and treeview.</p>

<p>pretty powerful!</p>
',1,1,0,5,5,0,'2009-10-16 18:25:38.237000',1,0,5);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-10-18 09:22:53.980000',0,'boy,george',20,0,0,'Do You Have To Be A Guy To Dress Up As Boy George',0,0,'<p>any ideas im a girl</p>
',20,0,'<p>any ideas im a girl</p>
',0,3,0,20,2,0,'2009-10-18 22:23:04.070000',1,7,15);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-10-18 09:23:34.373000',0,'boy,george',21,0,0,'Do You Have To Be A Guy To Dress Up As Boy George',0,0,'<p>any ideas im a girl</p>
',21,0,'<p>any ideas im a girl</p>
',0,3,0,21,2,0,'2009-10-18 22:25:02.687000',1,7,15);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-10-23 23:42:24.427000',5,'geneid,accession,mapping,conversion',22,493,13,'Gene Id Conversion Tool',0,4,'<p>Hey,</p>

<p>I was using DAVID (<a href=''http://david.abcc.ncifcrf.gov/conversion.jsp''>http://david.abcc.ncifcrf.gov/conversion.jsp</a>) to do the gene ID conversion, e.g.conversion between Agilent ID, Genebank accession id and Entrez gene ID, but I found the DAVID database is not updated. Does anyone know a better updataed conversion tool to do this job? Thanks! </p>
',22,0,'<p>Hey,</p>

<p>I was using DAVID (<a rel="nofollow" href="http://david.abcc.ncifcrf.gov/conversion.jsp">http://david.abcc.ncifcrf.gov/conversion.jsp</a>) to do the gene ID conversion, e.g.conversion between Agilent ID, Genebank accession id and Entrez gene ID, but I found the DAVID database is not updated. Does anyone know a better updataed conversion tool to do this job? Thanks! </p>
',0,1,2,22,59,0,'2010-03-03 20:03:06.160000',1,42,16);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-10-24 01:46:45.407000',0,'',23,0,0,'A: Gene Id Conversion Tool',0,0,'<p>I don''t know of a direct solution myself, but this is a topic that may be of interest for the biological data analysis class that I am teaching. </p>

<p>If you specify the organism/genomic builds that you are interested in we may be able to generate a full translation list as an in class example or a homework. I was planning on covering an <code>Affymetrix ID</code> to <code>Genebank example</code> anyhow.</p>
',22,0,'<p>I don''t know of a direct solution myself, but this is a topic that may be of interest for the biological data analysis class that I am teaching. </p>

<p>If you specify the organism/genomic builds that you are interested in we may be able to generate a full translation list as an in class example or a homework. I was planning on covering an <code>Affymetrix ID</code> to <code>Genebank example</code> anyhow.</p>
',1,1,0,22,2,0,'2009-10-24 01:46:45.407000',1,0,2);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-12-01 13:13:53.637000',2,'shrimp,sequencing,short,aligner',24,32,4,'How To Set Shrimp Parameters For Best Sensitivity With 35Bp Colorspace Data?',0,1,'<p>Hi,</p>

<p>I have 35bp Solid colorspace sequencing data, and the actual sequences to be mapped are 20-25bp after removing the linker sequence.</p>

<p>I hope to find all the hits allowing no more than n mismatches (say n=3), not only the best hit.</p>

<p>I know there is a -M option to specify -M sensitivity,35bp. I wonder whether this setting will guarantee the best sensitivity in this case. Since my reads are only 20-25bp long, should I changed the default 4 spaced seeds to 3?</p>

<p>I''m new to SHRiMP, so I''d like to hear some suggestions on setting the parameters of SHRiMP to achieve the best sensitivity.</p>

<p>Thank you!</p>
',24,0,'<p>Hi,</p>

<p>I have 35bp Solid colorspace sequencing data, and the actual sequences to be mapped are 20-25bp after removing the linker sequence.</p>

<p>I hope to find all the hits allowing no more than n mismatches (say n=3), not only the best hit.</p>

<p>I know there is a -M option to specify -M sensitivity,35bp. I wonder whether this setting will guarantee the best sensitivity in this case. Since my reads are only 20-25bp long, should I changed the default 4 spaced seeds to 3?</p>

<p>I''m new to SHRiMP, so I''d like to hear some suggestions on setting the parameters of SHRiMP to achieve the best sensitivity.</p>

<p>Thank you!</p>
',0,1,0,24,14,0,'2009-12-01 13:13:53.637000',1,21,14);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-12-01 20:57:35.300000',0,'',25,0,0,'A: How To Set Shrimp Parameters For Best Sensitivity With 35Bp Colorspace Data?',0,2,'<p>I just read the SHRiMP manual again, but I think that their explanation about -M option may not be enough to answer your question. I usually use the "seed" mode by using -s, -n, and -w and the option -M is a new feature of the version 1.3.1, which I have never tried before.</p>

<p>I recommend for you to use the "seed" mode--the default would be good, but please adjust the -s option if you want more sensitivity. Always fast speed compensates sensitivity and the -M option seems to exist for this purpose.</p>

<p>Hope my message to be helpful for your project.</p>
',24,0,'<p>I just read the SHRiMP manual again, but I think that their explanation about -M option may not be enough to answer your question. I usually use the "seed" mode by using -s, -n, and -w and the option -M is a new feature of the version 1.3.1, which I have never tried before.</p>

<p>I recommend for you to use the "seed" mode--the default would be good, but please adjust the -s option if you want more sensitivity. Always fast speed compensates sensitivity and the -M option seems to exist for this purpose.</p>

<p>Hope my message to be helpful for your project.</p>
',1,1,0,24,18,0,'2009-12-01 20:57:35.300000',1,0,18);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-12-02 00:00:53.443000',0,'',26,0,0,'A: How To Set Shrimp Parameters For Best Sensitivity With 35Bp Colorspace Data?',0,1,'<blockquote>
  <p>Since my reads are only 20-25bp long,
  should I changed the default 4 spaced
  seeds to 3?</p>
</blockquote>

<p>while the shrimp manual says:</p>

<ul>
<li>We recommend using the default 4 seeds of weight 12 in most cases.</li>
</ul>

<p>you could try running on a smaller sample and see what happens. </p>
',24,0,'<blockquote>
  <p>Since my reads are only 20-25bp long,
  should I changed the default 4 spaced
  seeds to 3?</p>
</blockquote>

<p>while the shrimp manual says:</p>

<ul>
<li>We recommend using the default 4 seeds of weight 12 in most cases.</li>
</ul>

<p>you could try running on a smaller sample and see what happens. </p>
',1,1,0,24,2,0,'2009-12-02 00:00:53.443000',1,0,2);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2009-12-09 03:45:54.547000',0,'',27,0,0,'A: Gene Id Conversion Tool',0,2,'<p>The following link has a list of ID conversion tools:</p>

<p><a href=''http://hum-molgen.org/NewsGen/08-2009/000020.html''>http://hum-molgen.org/NewsGen/08-2009/000020.html</a></p>
',22,0,'<p>The following link has a list of ID conversion tools:</p>

<p><a rel="nofollow" href="http://hum-molgen.org/NewsGen/08-2009/000020.html">http://hum-molgen.org/NewsGen/08-2009/000020.html</a></p>
',1,1,0,22,19,0,'2009-12-09 03:45:54.547000',1,0,19);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-01-13 18:59:22.603000',1,'meme,sge,motif,motif,compilation',28,59,1,'Tips On Compiling And Using Meme 4.3 With A Sun Grid Engine Computation Cluster',0,1,'<p>Has anyone compiled and used MEME 4.x for use in a parallel computation environment, based upon operation with a Sun Grid Engine (SGE) cluster?</p>

<p>I can compile the suite and its tests pass. However, when I attempt to use the <code>-p n</code> option, to specify <code>n</code> computation nodes, I get several error messages:</p>

<pre><code>/gridware/codine/util/arch: Command not found.
/gridware/codine/util/arch: Command not found.
/gridware/codine/util/arch: Command not found.
/gridware/codine/util/arch: Command not found.
1: Command not found.
</code></pre>

<p>We do not have <code>/gridware/codine/util/arch</code>, but we do have <code>/gridengine/sgi/util/arch</code>.</p>

<p>I tried looking around MEME''s source code, particularly at <code>meme.c</code> and <code>mp.h</code>, but there are no references to these paths.</p>

<p>I''m wondering if I am missing makefile directives. Here is my <code>./configure</code> statement:</p>

<pre><code>./configure --prefix=/home/areynolds/proj/meme/meme_4.3.0_build --with-url="http://meme.nbcr.net/meme" --enable-openmp --enable-debug
</code></pre>

<p>Is MPI a requirement; are there directives I am missing for MPI?</p>

<p>Thank you for any advice.</p>
',28,0,'<p>Has anyone compiled and used MEME 4.x for use in a parallel computation environment, based upon operation with a Sun Grid Engine (SGE) cluster?</p>

<p>I can compile the suite and its tests pass. However, when I attempt to use the <code>-p n</code> option, to specify <code>n</code> computation nodes, I get several error messages:</p>

<pre><code>/gridware/codine/util/arch: Command not found.
/gridware/codine/util/arch: Command not found.
/gridware/codine/util/arch: Command not found.
/gridware/codine/util/arch: Command not found.
1: Command not found.
</code></pre>

<p>We do not have <code>/gridware/codine/util/arch</code>, but we do have <code>/gridengine/sgi/util/arch</code>.</p>

<p>I tried looking around MEME''s source code, particularly at <code>meme.c</code> and <code>mp.h</code>, but there are no references to these paths.</p>

<p>I''m wondering if I am missing makefile directives. Here is my <code>./configure</code> statement:</p>

<pre><code>./configure --prefix=/home/areynolds/proj/meme/meme_4.3.0_build --with-url="<a rel="nofollow" href="http://meme.nbcr.net/meme">http://meme.nbcr.net/meme</a>" --enable-openmp --enable-debug
</code></pre>

<p>Is MPI a requirement; are there directives I am missing for MPI?</p>

<p>Thank you for any advice.</p>
',0,1,0,28,20,0,'2010-01-13 18:59:22.603000',1,14,20);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-01-14 03:17:50.023000',0,'',29,0,0,'A: Tips On Compiling And Using Meme 4.3 With A Sun Grid Engine Computation Cluster',0,0,'<p>This may not be overly useful but it very much sounds like a configuration problem.</p>

<p>Usually there is  configure flag that needs to be set to point to the libraries, something like:</p>

<pre><code>--with-mpidir=MPIDIR
--with-mpicc=MPICC
</code></pre>

<p>It also appears that the MEME suite does not support Open MPI (as per <a href=''http://meme.sdsc.edu/meme4/meme-install.html''>install notes</a>).</p>

<p>I would also recommend posting on the MEME user forum:</p>

<p><a href=''https://www.nbcr.net/forum/viewforum.php?f=5''>https://www.nbcr.net/forum/viewforum.php?f=5</a></p>
',28,0,'<p>This may not be overly useful but it very much sounds like a configuration problem.</p>

<p>Usually there is  configure flag that needs to be set to point to the libraries, something like:</p>

<pre><code>--with-mpidir=MPIDIR
--with-mpicc=MPICC
</code></pre>

<p>It also appears that the MEME suite does not support Open MPI (as per <a rel="nofollow" href="http://meme.sdsc.edu/meme4/meme-install.html">install notes</a>).</p>

<p>I would also recommend posting on the MEME user forum:</p>

<p><a rel="nofollow" href="https://www.nbcr.net/forum/viewforum.php?f=5">https://www.nbcr.net/forum/viewforum.php?f=5</a></p>
',1,1,0,28,2,0,'2010-01-14 03:17:50.023000',1,0,2);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-01-15 14:48:14.660000',0,'',30,0,0,'A: How Do I Convert From Bed Format To Gff Format?',0,2,'<p>Here''s a Perl script I wrote if you wanted to do something local. </p>

<p>There''s some code in there for translating yeast chromosome names that can be removed, if not needed. I also used a <code>Site</code> feature in the GFF file as the region ID, which might also need tweaking, depending on what features you''re interested in.</p>

<pre><code>#!/usr/bin/perl -w

use strict;
use Bio::Tools::GFF;
use feature qw(say switch);

my $gffio = Bio::Tools::GFF-&gt;new(-fh =&gt; \*STDIN, -gff_version =&gt; 2);
my $feature;

while ($feature = $gffio-&gt;next_feature()) {
    # print $gffio-&gt;gff_string($feature)."\n";

    # cf. http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml
    my $seq_id = $feature-&gt;seq_id();   
    my $start = $feature-&gt;start() - 1;
    my $end = $feature-&gt;end();
    my $strand = $feature-&gt;strand();
    my @sites = $feature-&gt;get_tag_values(''Site'');

    # translate strand
    given ( $strand ) {
        when ($_ == 1)  { $strand = "+"; }
        when ($_ == -1) { $strand = "-"; }
    }

    # translate yeast chromosome to UCSC browser-readable chromosome
    # cf. http://www.yeastgenome.org/sgdpub/Saccharomyces_cerevisiae.pdf
    given ( $seq_id ) {
        when ( $_ eq "I" )    { $seq_id = "chr1"; }
        when ( $_ eq "II" )   { $seq_id = "chr2"; }
        when ( $_ eq "III" )  { $seq_id = "chr3"; }
        when ( $_ eq "IV" )   { $seq_id = "chr4"; }
        when ( $_ eq "V" )    { $seq_id = "chr5"; }
        when ( $_ eq "VI" )   { $seq_id = "chr6"; }
        when ( $_ eq "VII" )  { $seq_id = "chr7"; }
        when ( $_ eq "VIII" ) { $seq_id = "chr8"; }
        when ( $_ eq "IX" )   { $seq_id = "chr9"; }
        when ( $_ eq "X" )    { $seq_id = "chr10"; }
        when ( $_ eq "XI" )   { $seq_id = "chr11"; }
        when ( $_ eq "XII" )  { $seq_id = "chr12"; }
        when ( $_ eq "XIII" ) { $seq_id = "chr13"; }
        when ( $_ eq "XIV" )  { $seq_id = "chr14"; }
        when ( $_ eq "XV" )   { $seq_id = "chr15"; }
        when ( $_ eq "XVI" )  { $seq_id = "chr16"; }
        default { }
    }

    # output
    print "$seq_id\t$start\t$end\t$sites[0]\t0.0\t$strand\n";
}
$gffio-&gt;close();
</code></pre>

<p>To use it:</p>

<pre><code>gff2bed.pl &lt; data.gff &gt; data.bed
</code></pre>
',2,0,'<p>Here''s a Perl script I wrote if you wanted to do something local. </p>

<p>There''s some code in there for translating yeast chromosome names that can be removed, if not needed. I also used a <code>Site</code> feature in the GFF file as the region ID, which might also need tweaking, depending on what features you''re interested in.</p>

<pre><code>#!/usr/bin/perl -w

use strict;
use Bio::Tools::GFF;
use feature qw(say switch);

my $gffio = Bio::Tools::GFF-&gt;new(-fh =&gt; \*STDIN, -gff_version =&gt; 2);
my $feature;

while ($feature = $gffio-&gt;next_feature()) {
    # print $gffio-&gt;gff_string($feature)."\n";

    # cf. <a rel="nofollow" href="http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml">http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml</a>
    my $seq_id = $feature-&gt;seq_id();   
    my $start = $feature-&gt;start() - 1;
    my $end = $feature-&gt;end();
    my $strand = $feature-&gt;strand();
    my @sites = $feature-&gt;get_tag_values(''Site'');

    # translate strand
    given ( $strand ) {
        when ($_ == 1)  { $strand = "+"; }
        when ($_ == -1) { $strand = "-"; }
    }

    # translate yeast chromosome to UCSC browser-readable chromosome
    # cf. <a rel="nofollow" href="http://www.yeastgenome.org/sgdpub/Saccharomyces_cerevisiae.pdf">http://www.yeastgenome.org/sgdpub/Saccharomyces_cerevisiae.pdf</a>
    given ( $seq_id ) {
        when ( $_ eq "I" )    { $seq_id = "chr1"; }
        when ( $_ eq "II" )   { $seq_id = "chr2"; }
        when ( $_ eq "III" )  { $seq_id = "chr3"; }
        when ( $_ eq "IV" )   { $seq_id = "chr4"; }
        when ( $_ eq "V" )    { $seq_id = "chr5"; }
        when ( $_ eq "VI" )   { $seq_id = "chr6"; }
        when ( $_ eq "VII" )  { $seq_id = "chr7"; }
        when ( $_ eq "VIII" ) { $seq_id = "chr8"; }
        when ( $_ eq "IX" )   { $seq_id = "chr9"; }
        when ( $_ eq "X" )    { $seq_id = "chr10"; }
        when ( $_ eq "XI" )   { $seq_id = "chr11"; }
        when ( $_ eq "XII" )  { $seq_id = "chr12"; }
        when ( $_ eq "XIII" ) { $seq_id = "chr13"; }
        when ( $_ eq "XIV" )  { $seq_id = "chr14"; }
        when ( $_ eq "XV" )   { $seq_id = "chr15"; }
        when ( $_ eq "XVI" )  { $seq_id = "chr16"; }
        default { }
    }

    # output
    print "$seq_id\t$start\t$end\t$sites[0]\t0.0\t$strand\n";
}
$gffio-&gt;close();
</code></pre>

<p>To use it:</p>

<pre><code>gff2bed.pl &lt; data.gff &gt; data.bed
</code></pre>
',1,1,0,2,20,0,'2010-01-15 14:48:14.660000',1,0,20);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-01-22 09:14:17.380000',2,'solid,rna,galaxy',31,79,4,'How Do I Map, Align, And Plot My Solid Results?',0,2,'<p>Hi, I recently performed an RNA immunoprecipitation followed by SOLiD sequencing (50 bp fragmented reads). I haven''t received my first SOLiD sequencing results yet, but I was told I should have them soon. I''ve tried doing my own research on how to map, align, and plot my results but I don''t have a concrete workflow as to how I will analyze my results yet. I have very little experience doing any programming and would prefer to use galaxy. There are labs on my campus I can go to to get my color space data mapped, but I would like to do things myself. Is there a way on galaxy (or another program) to convert my color space data to sequence, then map those reads to the yeast transcriptome and analyze it? Even if you can''t answer my question directly I''d appreciate any tips from anyone who has worked with RNA-seq data already.    </p>

<p>Thanks in advance </p>
',31,0,'<p>Hi, I recently performed an RNA immunoprecipitation followed by SOLiD sequencing (50 bp fragmented reads). I haven''t received my first SOLiD sequencing results yet, but I was told I should have them soon. I''ve tried doing my own research on how to map, align, and plot my results but I don''t have a concrete workflow as to how I will analyze my results yet. I have very little experience doing any programming and would prefer to use galaxy. There are labs on my campus I can go to to get my color space data mapped, but I would like to do things myself. Is there a way on galaxy (or another program) to convert my color space data to sequence, then map those reads to the yeast transcriptome and analyze it? Even if you can''t answer my question directly I''d appreciate any tips from anyone who has worked with RNA-seq data already.    </p>

<p>Thanks in advance </p>
',0,1,0,31,4,0,'2010-01-22 09:14:17.380000',1,17,4);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-01-22 21:13:42.790000',0,'',32,0,0,'A: How Do I Map, Align, And Plot My Solid Results?',0,1,'<p>Personally I would advise that if you know someone who can partially perform the task you should have them do it, and ask them to explain and show it to you how they''ve done it.</p>

<p>The task at hand is complex. The solution always depends immensely on the particulars of the problem, moreover you will be facing myriads of frustrating limitations, errors and problems.</p>

<p>Learning directly from someone who has done it, establishing a personal rapport with them will allow you to ease into this problem domain. In fact when you are finished mapping your RNA - your are still likely to be far from being done - yet you might have expanded a lot of energy and excitement. </p>
',31,0,'<p>Personally I would advise that if you know someone who can partially perform the task you should have them do it, and ask them to explain and show it to you how they''ve done it.</p>

<p>The task at hand is complex. The solution always depends immensely on the particulars of the problem, moreover you will be facing myriads of frustrating limitations, errors and problems.</p>

<p>Learning directly from someone who has done it, establishing a personal rapport with them will allow you to ease into this problem domain. In fact when you are finished mapping your RNA - your are still likely to be far from being done - yet you might have expanded a lot of energy and excitement. </p>
',1,1,0,31,2,0,'2010-01-22 21:13:42.790000',1,0,2);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-01-26 21:14:38.623000',4,'general,subjective,os',33,307,9,'Which Operating System Do You Prefer For Bioinformatics?',0,2,'<p>So, you will probably hate me for asking this question here, as there are lot of forum and blog posts on internet about it and it is also a very subjective question.</p>

<p>However, it may be a starting point for a good discussion, if we don''t flame... Which operating system do you usually use for your work? Did you install it by yourself, and do you have administrative rights on it, or is there any IT administrator in your lab? </p>
',33,0,'<p>So, you will probably hate me for asking this question here, as there are lot of forum and blog posts on internet about it and it is also a very subjective question.</p>

<p>However, it may be a starting point for a good discussion, if we don''t flame... Which operating system do you usually use for your work? Did you install it by yourself, and do you have administrative rights on it, or is there any IT administrator in your lab? </p>
',0,1,0,33,23,0,'2010-03-03 03:16:12.330000',1,28,23);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-01-26 21:45:06.737000',4,'subjective,programming,languages',34,486,14,'Which Are The Best Programming Languages To Study For A Bioinformatician?',0,3,'<p>This is also a very classic question: Which is your favorite programming language in bioinformatics? Which languages would you recommend to a student wishing to enter the world of bioinformatics?</p>

<p>This topic has already been discussed on the Internet, but I think it would be nice to discuss it here. Here there are some links to previous polls and discussions:</p>

<ul>
<li><a href=''http://www.bioinformatics.org/poll/index.php?dispid=17''>Bioinformatics.org poll</a></li>
<li><a href=''http://openwetware.org/wiki/Biogang:Projects/Bioinformatics_Career_Survey_2008''>Bioinformatics Career survey 2008 by Michael Barton</a></li>
</ul>
',34,0,'<p>This is also a very classic question: Which is your favorite programming language in bioinformatics? Which languages would you recommend to a student wishing to enter the world of bioinformatics?</p>

<p>This topic has already been discussed on the Internet, but I think it would be nice to discuss it here. Here there are some links to previous polls and discussions:</p>

<ul>
<li><a rel="nofollow" href="http://www.bioinformatics.org/poll/index.php?dispid=17">Bioinformatics.org poll</a></li>
<li><a rel="nofollow" href="http://openwetware.org/wiki/Biogang:Projects/Bioinformatics_Career_Survey_2008">Bioinformatics Career survey 2008 by Michael Barton</a></li>
</ul>
',0,1,0,34,23,0,'2010-04-15 16:08:30.167000',1,28,23);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-01-27 02:23:45.583000',0,'',35,0,0,'A: Which Operating System Do You Prefer For Bioinformatics?',0,1,'<p>Often people are limited to their choices by factors outside of their control. One lab that I work with requires the use of Mac computers another is using Windows mostly. Large scale computations seem to be best suited for Linux systems.</p>

<p>Luckily there is a migration towards unified capabilities across all platforms. Installing Cygwin on Windows allows us to tap into the power of Unix, while Linux distros have advanced graphical user interfaces like Windows and Macs.</p>

<p>From my own observations of non technical people, the installation of new and interdependent software packages seems to be the most difficult on Mac computers and easiest on Windows due to the computational architecture that makes all Windows computers identical. </p>
',33,0,'<p>Often people are limited to their choices by factors outside of their control. One lab that I work with requires the use of Mac computers another is using Windows mostly. Large scale computations seem to be best suited for Linux systems.</p>

<p>Luckily there is a migration towards unified capabilities across all platforms. Installing Cygwin on Windows allows us to tap into the power of Unix, while Linux distros have advanced graphical user interfaces like Windows and Macs.</p>

<p>From my own observations of non technical people, the installation of new and interdependent software packages seems to be the most difficult on Mac computers and easiest on Windows due to the computational architecture that makes all Windows computers identical. </p>
',1,1,0,33,2,0,'2010-01-27 02:23:45.583000',1,0,2);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-01-27 02:32:29.450000',0,'',36,0,0,'A: Which Are The Best Programming Languages To Study For A Bioinformatician?',0,1,'<p>It is important to be considerate and not characterize one particular approach negatively. My favorite quote is:</p>

<p><strong>Programming is pure thought.</strong></p>

<p>Hopefully everyone is able to pick an approach that matches their individual way of thinking. While I myself do not program in Perl, I consider it to be one of the most popular and powerful platforms for doing bioinformatics analysis. </p>
',34,0,'<p>It is important to be considerate and not characterize one particular approach negatively. My favorite quote is:</p>

<p><strong>Programming is pure thought.</strong></p>

<p>Hopefully everyone is able to pick an approach that matches their individual way of thinking. While I myself do not program in Perl, I consider it to be one of the most popular and powerful platforms for doing bioinformatics analysis. </p>
',1,1,0,34,2,0,'2010-01-27 02:32:29.450000',1,0,2);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-01-27 02:38:24.477000',0,'',37,0,0,'A: Which Operating System Do You Prefer For Bioinformatics?',0,1,'<p>Tips for installing software om Max OS X:</p>

<ul>
<li>install the Apple developer tools called <strong>Xcode</strong> <a href=''http://developer.apple.com/tools/xcode/''>http://developer.apple.com/tools/xcode/</a></li>
<li>install <strong>MacPorts</strong> from <a href=''http://www.macports.org/''>http://www.macports.org/</a></li>
</ul>

<p>You can now easily install everything from command line using the <code>port</code> command. List all available software</p>

<pre><code>port list
</code></pre>

<p>Install libraries and software. etc:</p>

<pre><code>port install &lt;some-library&gt;
</code></pre>
',33,0,'<p>Tips for installing software om Max OS X:</p>

<ul>
<li>install the Apple developer tools called <strong>Xcode</strong> <a rel="nofollow" href="http://developer.apple.com/tools/xcode/">http://developer.apple.com/tools/xcode/</a></li>
<li>install <strong>MacPorts</strong> from <a rel="nofollow" href="http://www.macports.org/">http://www.macports.org/</a></li>
</ul>

<p>You can now easily install everything from command line using the <code>port</code> command. List all available software</p>

<pre><code>port list
</code></pre>

<p>Install libraries and software. etc:</p>

<pre><code>port install &lt;some-library&gt;
</code></pre>
',1,1,0,33,2,0,'2010-01-27 02:38:24.477000',1,0,2);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-01-27 05:42:14.167000',0,'',38,0,0,'A: Which Are The Best Programming Languages To Study For A Bioinformatician?',0,1,'<p>Any programming language is good as long you know what you''re doing.</p>
',34,0,'<p>Any programming language is good as long you know what you''re doing.</p>
',1,1,0,34,25,0,'2010-01-27 05:42:14.167000',1,0,25);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-01-27 16:08:52.870000',0,'',39,0,0,'A: Gene Id Conversion Tool',0,2,'<p>You can also do it with the following services:</p>

<ul>
<li><a href=''http://www.uniprot.org/?tab=mapping''>uniprot</a> - Click on ''Id Mapping'' from the home page.</li>
<li><a href=''http://www.biomart.org/biomart/martview/''>biomart</a> - choose a database and a version, then put the ids you want to convert under Filters->Id List limit (select the proper input id in the menu), and then the output ids under ''Attributes''. Biomart is a general tool that enables you to extract a lot of different informations from databases - sequences, ontologies, transcripts, homologues - but maybe for converting gene ids is a bit too complex.</li>
<li><a href=''http://main.g2.bx.psu.edu/''>galaxy</a> - I can''t help too much about this here but I am sure it has a function for doing that - and many other things.</li>
</ul>
',22,0,'<p>You can also do it with the following services:</p>

<ul>
<li><a rel="nofollow" href="http://www.uniprot.org/?tab=mapping">uniprot</a> - Click on ''Id Mapping'' from the home page.</li>
<li><a rel="nofollow" href="http://www.biomart.org/biomart/martview/">biomart</a> - choose a database and a version, then put the ids you want to convert under Filters-&gt;Id List limit (select the proper input id in the menu), and then the output ids under ''Attributes''. Biomart is a general tool that enables you to extract a lot of different informations from databases - sequences, ontologies, transcripts, homologues - but maybe for converting gene ids is a bit too complex.</li>
<li><a rel="nofollow" href="http://main.g2.bx.psu.edu/">galaxy</a> - I can''t help too much about this here but I am sure it has a function for doing that - and many other things.</li>
</ul>
',1,1,0,22,23,0,'2010-01-27 16:08:52.870000',1,0,23);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-01-28 21:31:50.470000',0,'',40,0,0,'A: Finding Common Motifs In Sequences',0,1,'<p>Meme has been the first program to be published for doing that.
As an alternative you can find one of the <a href=''http://www.be.embnet.org/embosshelp/''>EMBOSS tools</a>; if you are scared by a terminal and want to do it from a web-based interface, you can use the EMBOSS tools from <a href=''http://main.g2.bx.psu.edu/''>galaxy</a></p>
',4,0,'<p>Meme has been the first program to be published for doing that.
As an alternative you can find one of the <a rel="nofollow" href="http://www.be.embnet.org/embosshelp/">EMBOSS tools</a>; if you are scared by a terminal and want to do it from a web-based interface, you can use the EMBOSS tools from <a rel="nofollow" href="http://main.g2.bx.psu.edu/">galaxy</a></p>
',1,1,0,4,23,0,'2010-01-28 21:31:50.470000',1,0,23);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-01-28 22:17:37.437000',2,'subjective,gene,go',41,119,6,'How Much Do You Trust Geneontology?',0,3,'<p><a href=''http://www.geneontology.org/''>GeneOntology</a> is a nice project to provide a standard terminology for genes and gene functions, to help avoid the use of synonyms and wrong spelling when describing a gene.</p>

<p>I have been using the GeneOntology for a while, but honestly I think that it contains many errors and that many terms have not enough terms associated. Moreover, the terminology they use is not always clear and there are some duplications.</p>

<p>It is frequent to read in article or in slideshows charts were the GO classification is used to infer the properties of a set of genes... But I wonder if the authors check the GO annotations that they use.</p>

<p>What is your experience about <a href=''http://www.geneontology.org/''>GO</a>?</p>
',41,0,'<p><a rel="nofollow" href="http://www.geneontology.org/">GeneOntology</a> is a nice project to provide a standard terminology for genes and gene functions, to help avoid the use of synonyms and wrong spelling when describing a gene.</p>

<p>I have been using the GeneOntology for a while, but honestly I think that it contains many errors and that many terms have not enough terms associated. Moreover, the terminology they use is not always clear and there are some duplications.</p>

<p>It is frequent to read in article or in slideshows charts were the GO classification is used to infer the properties of a set of genes... But I wonder if the authors check the GO annotations that they use.</p>

<p>What is your experience about <a rel="nofollow" href="http://www.geneontology.org/">GO</a>?</p>
',0,1,0,41,2,0,'2010-01-29 20:08:07.753000',1,21,23);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-01-28 22:22:04.203000',0,'',42,0,0,'A: Site Use Guidelines',0,0,'<p>The StackExchange websites have been designed for making questions related to programming and technical issues.</p>

<p>For example, for this reason, if you try to write a question which starts with ''What is your favorite experience...'' you get a disclaimer saying that ''your question seems to be probably subjective and it is likely to be closed''.</p>

<p>However, I think that it is very useful to make subjective and opinion-based questions on bioinformatics, as there are few places to do so... So, what is your policy? Will you accept subjective questions?</p>
',1,0,'<p>The StackExchange websites have been designed for making questions related to programming and technical issues.</p>

<p>For example, for this reason, if you try to write a question which starts with ''What is your favorite experience...'' you get a disclaimer saying that ''your question seems to be probably subjective and it is likely to be closed''.</p>

<p>However, I think that it is very useful to make subjective and opinion-based questions on bioinformatics, as there are few places to do so... So, what is your policy? Will you accept subjective questions?</p>
',1,1,0,1,23,0,'2010-01-28 22:22:04.203000',1,0,23);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-01-28 23:58:20.500000',0,'',43,0,0,'A: Which Are The Best Programming Languages To Study For A Bioinformatician?',0,7,'<p>The choice of a programming language is purely subjective, but when a student asks you which programming language he should start with, you have to make an answer, or at least provide some informations.</p>

<p>I think that a bioinformatician who studies <strong>R</strong> and at least two or three libraries (lattice/ggplot2, plyr) early can have an advantage, because he will be able to represent his data properly and obtain good results without too much effort. If your supervisor is not a computer scientist, he will be a lot more impressed by plots and charts than by programs, even if they are well written, with unittests etc.</p>

<p><strong>Python</strong> is a good programming language to learn as a general purpose tool. Its bigger advantages are its easy to read syntax, and its paradigm ''there is only one way to do it'', so the number of language keywords is reduced to the minimum, and two programs with the same function written by different people will be very similar (which is what doesn''t happen with perl). 
The negative points of python are that its CSV files reading/plotting interface is not ready yet (the best is pylab), so you must rely on R to produce nice plots.</p>

<p>Honestly I don''t like <strong>perl</strong>, because I think it can induce to many bad-behaviours in novel programmers. For example, in perl there are many similar constructs to accomplish the same objective: so, it is very difficult to understand a program written by someone else, because you have to known all the possible constructs and hope there are enough comments. It is already very difficult to reproduce a bioinformatician experiment, if you write your code in a difficult language it is a lot worst. 
Moreover, I know of many people who have been using perl for years, but that don''t even use functions, because it looks too complicated. How can it be? It looks very inefficient. 
The only good point of perl is its repositories, bioperl and CPAN; however, I know of people using perl that don''t even know of the existence of these, so I don''t understand why they keep going with perl.</p>

<p>Apart from programming language, is it very useful to learn the basic usage of <strong>gnu-make</strong>, or of a derivate. This program is very useful when you have lot of different scripts, as it allows you to define a pipeline in order to run them. 
Some basic <strong>bash commands</strong> may also be very useful if you work with a lot of flat files (head, sed, gawk, grep, ...)</p>
',34,0,'<p>The choice of a programming language is purely subjective, but when a student asks you which programming language he should start with, you have to make an answer, or at least provide some informations.</p>

<p>I think that a bioinformatician who studies <strong>R</strong> and at least two or three libraries (lattice/ggplot2, plyr) early can have an advantage, because he will be able to represent his data properly and obtain good results without too much effort. If your supervisor is not a computer scientist, he will be a lot more impressed by plots and charts than by programs, even if they are well written, with unittests etc.</p>

<p><strong>Python</strong> is a good programming language to learn as a general purpose tool. Its bigger advantages are its easy to read syntax, and its paradigm ''there is only one way to do it'', so the number of language keywords is reduced to the minimum, and two programs with the same function written by different people will be very similar (which is what doesn''t happen with perl). 
The negative points of python are that its CSV files reading/plotting interface is not ready yet (the best is pylab), so you must rely on R to produce nice plots.</p>

<p>Honestly I don''t like <strong>perl</strong>, because I think it can induce to many bad-behaviours in novel programmers. For example, in perl there are many similar constructs to accomplish the same objective: so, it is very difficult to understand a program written by someone else, because you have to known all the possible constructs and hope there are enough comments. It is already very difficult to reproduce a bioinformatician experiment, if you write your code in a difficult language it is a lot worst. 
Moreover, I know of many people who have been using perl for years, but that don''t even use functions, because it looks too complicated. How can it be? It looks very inefficient. 
The only good point of perl is its repositories, bioperl and CPAN; however, I know of people using perl that don''t even know of the existence of these, so I don''t understand why they keep going with perl.</p>

<p>Apart from programming language, is it very useful to learn the basic usage of <strong>gnu-make</strong>, or of a derivate. This program is very useful when you have lot of different scripts, as it allows you to define a pipeline in order to run them. 
Some basic <strong>bash commands</strong> may also be very useful if you work with a lot of flat files (head, sed, gawk, grep, ...)</p>
',1,1,0,34,23,0,'2010-01-28 23:58:20.500000',1,0,23);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-01-29 04:41:29.170000',0,'',44,0,0,'A: How Much Do You Trust Geneontology?',0,1,'<p>The GO terms and classifications are primarily an based on opinions and a human interpretation of a small group of people of what the current state of the knowledge is.Thus  are more subjective than say experimental measurements would be. </p>

<p>In fact it is surprising that it works at all; and it does indeed.  We just need to becareful not too read to much into it.</p>
',41,0,'<p>The GO terms and classifications are primarily an based on opinions and a human interpretation of a small group of people of what the current state of the knowledge is.Thus  are more subjective than say experimental measurements would be. </p>

<p>In fact it is surprising that it works at all; and it does indeed.  We just need to becareful not too read to much into it.</p>
',1,1,0,41,2,0,'2010-01-29 04:41:29.170000',1,0,2);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-01-29 11:50:26.020000',0,'',45,0,0,'A: How Much Do You Trust Geneontology?',0,2,'<p>In my experience it''s case by case. In other words just because you are getting significant p-values, does not mean the results are biologically significant. I once submitted clusters of microarray data and received a bunch of hits that were significant by p-value, but really didn''t have a theme. The GO terms I saw were from many different processes without an overall term (besides biological process) which linked them together. When I''ve looked at published GO terms searches I generally see a strong theme among many of the terms (however that doesn''t necessarily mean it has biological significance until tested empirically). So seeing themes among your terms may suggest higher significance, but it should make biological sense too. </p>
',41,0,'<p>In my experience it''s case by case. In other words just because you are getting significant p-values, does not mean the results are biologically significant. I once submitted clusters of microarray data and received a bunch of hits that were significant by p-value, but really didn''t have a theme. The GO terms I saw were from many different processes without an overall term (besides biological process) which linked them together. When I''ve looked at published GO terms searches I generally see a strong theme among many of the terms (however that doesn''t necessarily mean it has biological significance until tested empirically). So seeing themes among your terms may suggest higher significance, but it should make biological sense too. </p>
',1,1,0,41,4,0,'2010-01-29 11:50:26.020000',1,0,4);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-01-29 17:42:23.043000',2,'subjective,string,protein,interacti,ppi,pin',46,93,7,'What Is Your Experience With The String (Interactions) Database?',0,3,'<p>STRING is a database of predicted protein-protein interactions at EMBL. It cluster the results from many sources of protein-protein interactions databases, like Mint, etc.., and it also use the informations from KEGG-pathways and reactome, to provide the best annotations for the interactions of a protein.</p>

<p>I am a bit confused from the results that I see there, because when I look at the genes in the pathway I am studying, I see many errors and annotations that I don''t understand.</p>

<p>What is your experience with STRING? If you want to do me a favor, go there and try to see the interactions annotated for a gene that you know already. Do you see anything weird?</p>
',46,0,'<p>STRING is a database of predicted protein-protein interactions at EMBL. It cluster the results from many sources of protein-protein interactions databases, like Mint, etc.., and it also use the informations from KEGG-pathways and reactome, to provide the best annotations for the interactions of a protein.</p>

<p>I am a bit confused from the results that I see there, because when I look at the genes in the pathway I am studying, I see many errors and annotations that I don''t understand.</p>

<p>What is your experience with STRING? If you want to do me a favor, go there and try to see the interactions annotated for a gene that you know already. Do you see anything weird?</p>
',0,1,0,46,87,0,'2010-03-06 04:31:43.937000',1,21,23);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-01-29 20:06:06.180000',0,'',47,0,0,'A: What Is Your Experience With The String (Interactions) Database?',0,1,'<p>I have not used STRING in particular but I have worked with protein interactions before (DIP dataset). I recall that even experimentally produced protein-protein interactions may have very large false positive ratios  (as for false negatives, who knows?) Some papers claim that up to 50% of the interactions were spurious; and repeated experiments showed very small overlaps. Predictions may be even less reliable.</p>

<p>At the same time the DIP dataset performed substantially better if we only considered the interactions for which there were multiple sources of evidence, so that may be a strategy to consider in your case as well.</p>
',46,0,'<p>I have not used STRING in particular but I have worked with protein interactions before (DIP dataset). I recall that even experimentally produced protein-protein interactions may have very large false positive ratios  (as for false negatives, who knows?) Some papers claim that up to 50% of the interactions were spurious; and repeated experiments showed very small overlaps. Predictions may be even less reliable.</p>

<p>At the same time the DIP dataset performed substantially better if we only considered the interactions for which there were multiple sources of evidence, so that may be a strategy to consider in your case as well.</p>
',1,1,0,46,2,0,'2010-01-29 20:06:06.180000',1,0,2);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-12 23:33:06.820000',3,'sequence,protein,structure',48,157,10,'Where Can I Get The Secondary Structure Of A Protein?',0,2,'<p>As in the title... I have a protein and I would like to know its secundary structure.
I couldn''t find it in uniprot, althought I tought they had annotations for it there.
In the end I have used a predictor (<a href=''http://www.compbio.dundee.ac.uk/www-jpred''>jpred</a>) but there it should be a database somewhere.</p>
',48,0,'<p>As in the title... I have a protein and I would like to know its secundary structure.
I couldn''t find it in uniprot, althought I tought they had annotations for it there.
In the end I have used a predictor (<a rel="nofollow" href="http://www.compbio.dundee.ac.uk/www-jpred">jpred</a>) but there it should be a database somewhere.</p>
',0,1,0,48,23,0,'2010-02-27 15:58:45.557000',1,28,23);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-13 02:01:49.187000',0,'',49,0,0,'A: Where Can I Get The Secondary Structure Of A Protein?',0,1,'<p>Protein structure prediction is a complex issue that is likely to require multiple approaches. There are many methods/tools listed at the </p>

<ul>
<li><a href=''http://www.expasy.ch/''>Expert Protein Analysis System website</a></li>
</ul>
',48,0,'<p>Protein structure prediction is a complex issue that is likely to require multiple approaches. There are many methods/tools listed at the </p>

<ul>
<li><a rel="nofollow" href="http://www.expasy.ch/">Expert Protein Analysis System website</a></li>
</ul>
',1,1,0,48,2,0,'2010-02-13 02:01:49.187000',1,0,2);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-13 03:57:06.680000',0,'',50,0,0,'A: Where Can I Get The Secondary Structure Of A Protein?',0,3,'<p>I think you found the best answer yourself: use a predictor! There are several out there...</p>

<p>You suggest that there should be a Secondary Structure Database. I''m not sure that makes much sense, let me explain my point of view (which may not be that of everyone): most often, the data that is found in databases is the "state of knowledge" of the described object, based on experimentation.</p>

<p>That may be the case for secondary structures of proteins, but only in the case where the said proteins have been crystalized. In those cases, it is not only the secondary structures but also the tertiary structures (with the caveat that the crystal structure of a protein does not prove "all" states that a protein may take in real "dynamic" physiological conditions).</p>

<p>For all those proteins that have not been crystalized, then we can only rely on predictions. And I use them quite frequently: they are extremely useful! But as far as I know, no prediction is accepted as fact. They''re "educated guesses" that are often correct, but sometimes wrong. The results may differ from one prediction method to another. Also they change each time the algorithms are improved...</p>

<p>If there was a database of predicted secondary structures, people would likely take them for granted (make the equation prediction = fact) which would be quite "unscientific".</p>

<p>I think such a resource would be more of a hindrance than an asset to the scientific community...</p>
',48,0,'<p>I think you found the best answer yourself: use a predictor! There are several out there...</p>

<p>You suggest that there should be a Secondary Structure Database. I''m not sure that makes much sense, let me explain my point of view (which may not be that of everyone): most often, the data that is found in databases is the "state of knowledge" of the described object, based on experimentation.</p>

<p>That may be the case for secondary structures of proteins, but only in the case where the said proteins have been crystalized. In those cases, it is not only the secondary structures but also the tertiary structures (with the caveat that the crystal structure of a protein does not prove "all" states that a protein may take in real "dynamic" physiological conditions).</p>

<p>For all those proteins that have not been crystalized, then we can only rely on predictions. And I use them quite frequently: they are extremely useful! But as far as I know, no prediction is accepted as fact. They''re "educated guesses" that are often correct, but sometimes wrong. The results may differ from one prediction method to another. Also they change each time the algorithms are improved...</p>

<p>If there was a database of predicted secondary structures, people would likely take them for granted (make the equation prediction = fact) which would be quite "unscientific".</p>

<p>I think such a resource would be more of a hindrance than an asset to the scientific community...</p>
',1,1,0,48,26,0,'2010-02-13 03:57:06.680000',1,0,26);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-16 07:13:06.543000',1,'sequence,blast',51,55,2,'Turn Off Blast Search On Reverse Complement Strand In Blastn',0,1,'<p>I have a quick question:
How can I turn off search on reverse complement strand of my query nucleotide sequence in blastn?</p>

<p>For example, I don''t want ''GUAAAGCCAAAUCUUCGGUUA'' to be a hit when I use ''UAACCGAAGAUUUGGCUUUAC'' as the query.</p>

<p>Maybe I missed it when I read the man page, but I really appreciate it if someone can point out the parameter I should use.</p>

<p>Thanks!</p>
',51,0,'<p>I have a quick question:
How can I turn off search on reverse complement strand of my query nucleotide sequence in blastn?</p>

<p>For example, I don''t want ''GUAAAGCCAAAUCUUCGGUUA'' to be a hit when I use ''UAACCGAAGAUUUGGCUUUAC'' as the query.</p>

<p>Maybe I missed it when I read the man page, but I really appreciate it if someone can point out the parameter I should use.</p>

<p>Thanks!</p>
',0,1,0,51,2,0,'2010-02-17 01:20:31.193000',1,14,14);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-16 08:31:37.060000',0,'',52,0,0,'A: Turn Off Blast Search On Reverse Complement Strand In Blastn',0,1,'<p>The -S flag can select the strands:</p>

<pre><code>-S  Query strands to search against database 
    (for blast[nx], and tblastx) 3 is both, 1 is top, 2 is bottom [Integer]
</code></pre>
',51,0,'<p>The -S flag can select the strands:</p>

<pre><code>-S  Query strands to search against database 
    (for blast[nx], and tblastx) 3 is both, 1 is top, 2 is bottom [Integer]
</code></pre>
',1,1,0,51,2,0,'2010-02-16 08:31:37.060000',1,0,2);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-19 13:29:18.733000',1,'solid',53,53,2,'How To Do Quality Trimming Of Solid Reads In Colour Space?',0,1,'<p>The reads returned from the Solid sequencing provider are littered with dots and some bases have a negative quality value. Does anyone know if there is a good method to extract high quality regions from the reads without distorting the reading of bases in colour space?</p>
',53,0,'<p>The reads returned from the Solid sequencing provider are littered with dots and some bases have a negative quality value. Does anyone know if there is a good method to extract high quality regions from the reads without distorting the reading of bases in colour space?</p>
',0,1,0,53,27,0,'2010-02-19 13:29:18.733000',1,14,27);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-19 13:31:24.643000',0,'',54,0,0,'A: How Do I Map, Align, And Plot My Solid Results?',0,1,'<p>You can try <a href=''http://bio-bwa.sourceforge.net/''>BWA</a> as well:
<a href=''http://maq.sourceforge.net/bwa-man.shtml''>http://maq.sourceforge.net/bwa-man.shtml</a></p>
',31,0,'<p>You can try <a rel="nofollow" href="http://bio-bwa.sourceforge.net/">BWA</a> as well:
<a rel="nofollow" href="http://maq.sourceforge.net/bwa-man.shtml">http://maq.sourceforge.net/bwa-man.shtml</a></p>
',1,1,0,31,27,0,'2010-02-19 13:31:24.643000',1,0,27);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-19 19:33:45.533000',0,'',55,0,0,'A: How To Do Quality Trimming Of Solid Reads In Colour Space?',0,1,'<p>The <a href=''http://solidsoftwaretools.com/gf/project/saet/''>Solid Accuracy Enhancer Tool</a> might be useful for this.</p>
',53,0,'<p>The <a rel="nofollow" href="http://solidsoftwaretools.com/gf/project/saet/">Solid Accuracy Enhancer Tool</a> might be useful for this.</p>
',1,1,0,53,28,0,'2010-02-19 19:33:45.533000',1,0,28);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-21 22:13:39.320000',2,'sequence,ucsc,fasta',56,80,9,'How To Get The Sequence Of A Genomic Region From Ucsc?',0,4,'<p>Let''s say I want to download the fasta sequence of the region chr1:100000..200000 from the UCSC browser.
How do you do that? I can''t find a button to ''export to fasta'' in the UCSC genome browser. I think that the solution is to click on one of the tracks displayed, but I am not sure of which.
If I go to the Tables section, I can''t find a table with the fasta sequences among the many.</p>
',56,0,'<p>Let''s say I want to download the fasta sequence of the region chr1:100000..200000 from the UCSC browser.
How do you do that? I can''t find a button to ''export to fasta'' in the UCSC genome browser. I think that the solution is to click on one of the tracks displayed, but I am not sure of which.
If I go to the Tables section, I can''t find a table with the fasta sequences among the many.</p>
',0,1,0,56,2,0,'2010-02-22 21:40:46.360000',1,21,23);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-22 01:11:20.450000',0,'',57,0,0,'A: How To Get The Sequence Of A Genomic Region From Ucsc?',0,2,'<p>The Genome Browser is for visualization.</p>

<p>To get data in many formats use the <a href=''http://genome.ucsc.edu/cgi-bin/hgTables?org=human''>UCSC Table Browser</a> then select the output format of your choice.</p>

<p>You may also need to select the right <strong>group</strong> and <strong>track</strong> to get the data you want.</p>
',56,0,'<p>The Genome Browser is for visualization.</p>

<p>To get data in many formats use the <a rel="nofollow" href="http://genome.ucsc.edu/cgi-bin/hgTables?org=human">UCSC Table Browser</a> then select the output format of your choice.</p>

<p>You may also need to select the right <strong>group</strong> and <strong>track</strong> to get the data you want.</p>
',1,1,0,56,2,0,'2010-02-22 01:11:20.450000',1,0,2);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-25 21:39:15.467000',6,'general,subjective',58,205,12,'What Is The Best Way To Share Scripts Between Members Of A Lab?',0,1,'<p>One of the most awful problems in my group is avoiding to rewrite scripts that have been already written by others. Since we have different projects and we work with different data, everybody ends up writing its own scripts in his favorite programming language, and it is very frequent to waste an afternoon on writing a new program and then discover that your workmate already had a script to do that.</p>

<p>Apart from the most logical answer ("talk with your workmates"), we are thinking about having a common place to store our best scripts, and if possible work together on them.
It would be similar to an image library like this: <a href=''http://matplotlib.sourceforge.net/gallery.html''>http://matplotlib.sourceforge.net/gallery.html</a> , where to put the script and an example of its output (most of our scripts produce graphs), and if possible integrated with Git. </p>

<p>Do you have any idea? How to you cope with the problem in your lab?</p>
',58,0,'<p>One of the most awful problems in my group is avoiding to rewrite scripts that have been already written by others. Since we have different projects and we work with different data, everybody ends up writing its own scripts in his favorite programming language, and it is very frequent to waste an afternoon on writing a new program and then discover that your workmate already had a script to do that.</p>

<p>Apart from the most logical answer ("talk with your workmates"), we are thinking about having a common place to store our best scripts, and if possible work together on them.
It would be similar to an image library like this: <a rel="nofollow" href="http://matplotlib.sourceforge.net/gallery.html">http://matplotlib.sourceforge.net/gallery.html</a> , where to put the script and an example of its output (most of our scripts produce graphs), and if possible integrated with Git. </p>

<p>Do you have any idea? How to you cope with the problem in your lab?</p>
',0,1,0,58,23,0,'2010-02-25 21:39:15.467000',1,49,23);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-25 23:32:53.740000',0,'',59,0,0,'A: How To Get The Sequence Of A Genomic Region From Ucsc?',0,3,'<p>Use the <strong>DAS</strong> server:</p>

<p><a href=''http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr1:100000,200000''>http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr1:100000,200000</a></p>
',56,0,'<p>Use the <strong>DAS</strong> server:</p>

<p><a rel="nofollow" href="http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr1:100000,200000">http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr1:100000,200000</a></p>
',1,1,0,56,30,0,'2010-02-25 23:32:53.740000',1,0,30);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-25 23:51:28.290000',0,'',60,0,0,'A: Finding Common Motifs In Sequences',0,1,'<p>Some time ago I used SOMBRERO (<a href=''http://bioinf.nuigalway.ie/sombrero/download.html''>http://bioinf.nuigalway.ie/sombrero/download.html</a>) with a good degree of success on finding motifs in a very diverse set of sequences. They have a Mac version for download as well as parallel versions for Irix and Linux.</p>
',4,0,'<p>Some time ago I used SOMBRERO (<a rel="nofollow" href="http://bioinf.nuigalway.ie/sombrero/download.html">http://bioinf.nuigalway.ie/sombrero/download.html</a>) with a good degree of success on finding motifs in a very diverse set of sequences. They have a Mac version for download as well as parallel versions for Irix and Linux.</p>
',1,1,0,4,31,0,'2010-02-25 23:51:28.290000',1,0,31);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-25 23:54:49.080000',0,'',61,0,0,'A: What Is The Best Way To Share Scripts Between Members Of A Lab?',0,1,'<p>I would recommend you to setup a wiki for your group. If you do not have a server readily you can always use one of the many wiki services available for free like Wikispaces (www.wikispaces.com).</p>
',58,0,'<p>I would recommend you to setup a wiki for your group. If you do not have a server readily you can always use one of the many wiki services available for free like Wikispaces www.wikispaces.com).</p>
',1,1,0,58,31,0,'2010-02-25 23:54:49.080000',1,0,31);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-26 00:24:46.380000',0,'',62,0,0,'A: What Is The Best Way To Share Scripts Between Members Of A Lab?',0,2,'<p>Integrating with the source code management tool is essential, that way when code gets changed everyone can easily get the updated version. Wikis are also a good idea.</p>
',58,0,'<p>Integrating with the source code management tool is essential, that way when code gets changed everyone can easily get the updated version. Wikis are also a good idea.</p>
',1,1,0,58,2,0,'2010-02-26 00:24:46.380000',1,0,2);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-26 03:33:37.903000',0,'',63,0,0,'A: What Is Your Experience With The String (Interactions) Database?',0,3,'<p>I''ve been using STRING extensively, but not for protein-protein interactions work. STRING, as you note, is a bit of a mutt in terms of the different data sources it mines. Some that you''re missing include a broad literature-based search, as well as gene expression data sets. So if you''re interested primarily in physical interactions or any other single type of data source, STRING is a poor choice for your work. On the other hand, STRING does provide confidence scores for each association, as well as annotation for their data source types (with the license). So you can use those to filter out the interactions derived from data types you don''t want to see.</p>
',46,0,'<p>I''ve been using STRING extensively, but not for protein-protein interactions work. STRING, as you note, is a bit of a mutt in terms of the different data sources it mines. Some that you''re missing include a broad literature-based search, as well as gene expression data sets. So if you''re interested primarily in physical interactions or any other single type of data source, STRING is a poor choice for your work. On the other hand, STRING does provide confidence scores for each association, as well as annotation for their data source types (with the license). So you can use those to filter out the interactions derived from data types you don''t want to see.</p>
',1,1,0,46,35,0,'2010-02-26 03:33:37.903000',1,0,35);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-26 03:41:22.310000',0,'',64,0,0,'A: Which Are The Best Programming Languages To Study For A Bioinformatician?',0,2,'<p>Perl can be quite lovely if you choose to write it well. If you find yourself in need of writing some perl, I''d highly recommend getting the Perl Best Practices book and going through it to learn how to make your perl code not suck. Essential tools for helping with that are perlcritic and perltidy, both of which I have bound to quick keystrokes in my emacs cperl-mode so as to make sure my code is in reasonably good shape. There''s lots of blog articles out there about writing "Modern Perl" or "Enlightened Perl" that help make the language not just bearable but actually quite nice for a certain type of brain.</p>

<p>One thing that Perl does very well that no other language does is quick text processing on the command line. If you want to do some simple processing of a text file (which is pretty standard in this business), perl is a fantastic package to do so. Stringing together a set of UNIX utilities on a Linux system will usually have you running for a half dozen manpages looking for conflicting and unique switches, where with perl I find that there''s far less I have to remember to get the same effect. The book Minimal Perl goes in to this sort of thing in detail (perl as a better awk/sed/grep/etc) and I highly recommend having a look. At the very least, I''ve found that using perl in this fashion filled a hole in my toolkit that I didn''t even realize was there. R and Python can, of course, do this sort of thing too, but not nearly so well as Perl.</p>
',34,0,'<p>Perl can be quite lovely if you choose to write it well. If you find yourself in need of writing some perl, I''d highly recommend getting the Perl Best Practices book and going through it to learn how to make your perl code not suck. Essential tools for helping with that are perlcritic and perltidy, both of which I have bound to quick keystrokes in my emacs cperl-mode so as to make sure my code is in reasonably good shape. There''s lots of blog articles out there about writing "Modern Perl" or "Enlightened Perl" that help make the language not just bearable but actually quite nice for a certain type of brain.</p>

<p>One thing that Perl does very well that no other language does is quick text processing on the command line. If you want to do some simple processing of a text file (which is pretty standard in this business), perl is a fantastic package to do so. Stringing together a set of UNIX utilities on a Linux system will usually have you running for a half dozen manpages looking for conflicting and unique switches, where with perl I find that there''s far less I have to remember to get the same effect. The book Minimal Perl goes in to this sort of thing in detail (perl as a better awk/sed/grep/etc) and I highly recommend having a look. At the very least, I''ve found that using perl in this fashion filled a hole in my toolkit that I didn''t even realize was there. R and Python can, of course, do this sort of thing too, but not nearly so well as Perl.</p>
',1,1,0,34,35,0,'2010-02-26 03:41:22.310000',1,0,35);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-26 03:50:25.157000',0,'',65,0,0,'A: Which Operating System Do You Prefer For Bioinformatics?',0,2,'<p>My tip: install Cygwin if you are using Windows </p>
',33,0,'<p>My tip: install Cygwin if you are using Windows </p>
',1,1,0,33,33,0,'2010-02-26 03:50:25.157000',1,0,33);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-26 03:52:10.690000',0,'',66,0,0,'A: Site Use Guidelines',0,2,'<p>Who is running this site?</p>
',1,0,'<p>Who is running this site?</p>
',1,1,0,1,33,0,'2010-02-26 03:52:10.690000',1,0,33);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-26 15:51:33.997000',0,'',67,0,0,'A: Which Operating System Do You Prefer For Bioinformatics?',0,3,'<p>All of the 3 major platforms have their advantages, and I use all 3 practically every day. Mac OS X is my primary desktop OS, for a number of reasons, but mostly because I just seem more productive using it than any of the alternatives. All of my coding work is done over SSH on Linux (almost exclusively Ubuntu) servers. The power of Aptitude package management, and the robustness of this platform means that there really is no other choice for this kind of work. Finally I run Windows 7 on my netbook, because it is an excellent OS for that platform, and enables me to do everything I want that machine to be capable of, note-taking, blog writing, as a display machine for Powerpoint etc. It is also useful to have Internet Explorer kicking around somewhere for compatability testing.</p>

<p>I wouldn''t consider using any machine that I didn''t have admin rights on for work purposes, if I have to jump through hoops to get stuff installed, it just slows me down too much. This is another reason for using OS X for my primary desktop, it allows me to escape the University''s "Common Desktop" policy for Windows PCs, which would take control of my computer out of my hands.</p>
',33,0,'<p>All of the 3 major platforms have their advantages, and I use all 3 practically every day. Mac OS X is my primary desktop OS, for a number of reasons, but mostly because I just seem more productive using it than any of the alternatives. All of my coding work is done over SSH on Linux (almost exclusively Ubuntu) servers. The power of Aptitude package management, and the robustness of this platform means that there really is no other choice for this kind of work. Finally I run Windows 7 on my netbook, because it is an excellent OS for that platform, and enables me to do everything I want that machine to be capable of, note-taking, blog writing, as a display machine for Powerpoint etc. It is also useful to have Internet Explorer kicking around somewhere for compatability testing.</p>

<p>I wouldn''t consider using any machine that I didn''t have admin rights on for work purposes, if I have to jump through hoops to get stuff installed, it just slows me down too much. This is another reason for using OS X for my primary desktop, it allows me to escape the University''s "Common Desktop" policy for Windows PCs, which would take control of my computer out of my hands.</p>
',1,1,0,33,38,0,'2010-02-26 15:56:37.950000',1,0,38);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-26 16:07:14.700000',0,'',68,0,0,'A: What Is The Best Way To Share Scripts Between Members Of A Lab?',0,2,'<p>If you want to see the code, but also store associated information, such as expected outputs etc, then a wiki probably is the best choice (we prefer <a href=''http://www.dokuwiki.org/dokuwiki''>DokuWiki</a> here), although this would involve a lot of manual effort to document each script. </p>

<p>Use of a site such as GitHub would give you version control + a handy place to read code, although it is not free to host private repositories there, which I guess is what the majority of labs would require. </p>

<p>If privacy is not a concern, then I would consider GitHub <a href=''http://gist.github.com/''>gists</a> for code, which can then be embedded in a <a href=''http://posterous.com/''>Posterous</a> blog for comments. Posterous automatically unfolds Gist URLs into code samples in blog posts, so then you can annotate them easily. This would be a lot less manual effort than a wiki.</p>
',58,0,'<p>If you want to see the code, but also store associated information, such as expected outputs etc, then a wiki probably is the best choice (we prefer <a rel="nofollow" href="http://www.dokuwiki.org/dokuwiki">DokuWiki</a> here), although this would involve a lot of manual effort to document each script. </p>

<p>Use of a site such as GitHub would give you version control + a handy place to read code, although it is not free to host private repositories there, which I guess is what the majority of labs would require. </p>

<p>If privacy is not a concern, then I would consider GitHub <a rel="nofollow" href="http://gist.github.com/">gists</a> for code, which can then be embedded in a <a rel="nofollow" href="http://posterous.com/">Posterous</a> blog for comments. Posterous automatically unfolds Gist URLs into code samples in blog posts, so then you can annotate them easily. This would be a lot less manual effort than a wiki.</p>
',1,1,0,58,38,0,'2010-02-26 16:07:14.700000',1,0,38);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-26 18:50:17.577000',5,'hdf,biohdf,hdf,storage',69,314,20,'Using Hdf5 To Store  Bio-Data',0,8,'<p>Hi all,
has anobody ever used the <a href=''http://www.hdfgroup.org/HDF5/''>HDF5 API</a> to store some biological data (genotypes...). I know about this <a href=''http://www.geospiza.com/finchtalk/2008/03/genotyping-with-hdf.html''>kind of reference</a> (BioHDF...)  but I''m looking for some <strong>source code</strong> I could browse to understand how I can access data faster.</p>

<p>Pierre</p>

<p>PS: hum, I''m a new user. I''m not allowed to add the following tags: storage database hdf5 source code </p>
',69,0,'<p>Hi all,
has anobody ever used the <a rel="nofollow" href="http://www.hdfgroup.org/HDF5/">HDF5 API</a> to store some biological data (genotypes...). I know about this <a rel="nofollow" href="http://www.geospiza.com/finchtalk/2008/03/genotyping-with-hdf.html">kind of reference</a> (BioHDF...)  but I''m looking for some <strong>source code</strong> I could browse to understand how I can access data faster.</p>

<p>Pierre</p>

<p>PS: hum, I''m a new user. I''m not allowed to add the following tags: storage database hdf5 source code </p>
',0,1,0,69,30,0,'2010-03-03 19:38:58.783000',1,35,30);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-26 19:15:16.893000',0,'',70,0,0,'A: What Is The Best Way To Share Scripts Between Members Of A Lab?',0,2,'<p>You might also want to setup a simple snippets database. Navysnip application by Jason Strutz is easy to install and run if you have ruby and rubyonrails installed.</p>

<p>git clone git://github.com/navyrain/navysnip.git
  cd navysnip
  sudo rake gems:install
  rake db:migrate
  ruby script/server
Then visit your app at <a href=''http://localhost:3000''>http://localhost:3000</a></p>

<p>check out <a href=''http://github.com/navyrain/navysnip''>http://github.com/navyrain/navysnip</a>  for complete details</p>
',58,0,'<p>You might also want to setup a simple snippets database. Navysnip application by Jason Strutz is easy to install and run if you have ruby and rubyonrails installed.</p>

<p>git clone git://github.com/navyrain/navysnip.git
  cd navysnip
  sudo rake gems:install
  rake db:migrate
  ruby script/server
Then visit your app at <a rel="nofollow" href="http://localhost:3000">http://localhost:3000</a></p>

<p>check out <a rel="nofollow" href="http://github.com/navyrain/navysnip">http://github.com/navyrain/navysnip</a>  for complete details</p>
',1,1,0,58,40,0,'2010-02-26 19:15:16.893000',1,0,40);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-26 19:47:41.697000',0,'',71,0,0,'A: Using Hdf5 To Store  Bio-Data',0,2,'<p>Hello Pierre!</p>

<p>I have been talking with the BioHDF guys and from what they tell me, their work will be centered around a number of command-line APIs, written in C, that will address some areas of usage which for now do not seem to overlap. </p>

<p>I have seen this example on their site:
<a href=''http://www.hdfgroup.org/projects/biohdf/biohdf_tools.html''>http://www.hdfgroup.org/projects/biohdf/biohdf_tools.html</a>
Don''t know if that helps.</p>

<p>I have been talking with them to see if we can achieve an API for saving genotype data. Don''t know yet where that will lead me.</p>

<p>If you are looking for something more versatile, you will probably have to delve in the official HDF5 C code ( <a href=''http://www.hdfgroup.org/HDF5/Tutor/''>http://www.hdfgroup.org/HDF5/Tutor/</a> ), which seems to be the only one that offers all the functionality and goodies of that impressive storage system. </p>
',69,0,'<p>Hello Pierre!</p>

<p>I have been talking with the BioHDF guys and from what they tell me, their work will be centered around a number of command-line APIs, written in C, that will address some areas of usage which for now do not seem to overlap. </p>

<p>I have seen this example on their site:
<a rel="nofollow" href="http://www.hdfgroup.org/projects/biohdf/biohdf_tools.html">http://www.hdfgroup.org/projects/biohdf/biohdf_tools.html</a>
Don''t know if that helps.</p>

<p>I have been talking with them to see if we can achieve an API for saving genotype data. Don''t know yet where that will lead me.</p>

<p>If you are looking for something more versatile, you will probably have to delve in the official HDF5 C code ( <a rel="nofollow" href="http://www.hdfgroup.org/HDF5/Tutor/">http://www.hdfgroup.org/HDF5/Tutor/</a> ), which seems to be the only one that offers all the functionality and goodies of that impressive storage system. </p>
',1,1,0,69,42,0,'2010-02-26 19:47:41.697000',1,0,42);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-26 19:54:03.447000',0,'',72,0,0,'A: Using Hdf5 To Store  Bio-Data',0,1,'<p>Unfortunately I don''t have any example to shows you yet.
I don''t know how to program in C/C++ so I have been looking at two hdf5 wrappers in python, <a href=''http://www.pytables.org/moin''>PyTables</a> and <a href=''http://h5py.alfven.org/''>H5PY</a>.</p>

<p>PyTables has a database-like approach in which HDF5 is used as a sort of hierarchical database, in which a column can be a table itself, allowing to store nested data. For example, you have a table called ''SNPs'' with two columns, ''id'' and ''genotypes''; the column ''genotypes'' contains a nested table, with the columns ''individual'' and ''genotype''; and so on.</p>

<p>H5Py is basically a re-implementation of numpy''s arrays, so you can store and access arrays/matrixes as you would do with numpy (it is similar to arrays and matrixes in matlab, R, and any other language with this data type) and they are stored in an HDF5 file so the access is faster.</p>
',69,0,'<p>Unfortunately I don''t have any example to shows you yet.
I don''t know how to program in C/C++ so I have been looking at two hdf5 wrappers in python, <a rel="nofollow" href="http://www.pytables.org/moin">PyTables</a> and <a rel="nofollow" href="http://h5py.alfven.org/">H5PY</a>.</p>

<p>PyTables has a database-like approach in which HDF5 is used as a sort of hierarchical database, in which a column can be a table itself, allowing to store nested data. For example, you have a table called ''SNPs'' with two columns, ''id'' and ''genotypes''; the column ''genotypes'' contains a nested table, with the columns ''individual'' and ''genotype''; and so on.</p>

<p>H5Py is basically a re-implementation of numpy''s arrays, so you can store and access arrays/matrixes as you would do with numpy (it is similar to arrays and matrixes in matlab, R, and any other language with this data type) and they are stored in an HDF5 file so the access is faster.</p>
',1,1,0,69,23,0,'2010-02-26 19:54:03.447000',1,0,23);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-26 19:56:22.753000',0,'',73,0,0,'A: Using Hdf5 To Store  Bio-Data',0,3,'<p>In the <a href=''http://code.google.com/p/genetrack/''>GeneTrack</a> software we have used HDF to store values for each genomic base. Its main advantage over other storage systems was that it was able to return consecutive values with minimal overhead. </p>

<p>For example it is <em>extremely fast</em> (ms) in retrieving say 100,000 consecutive values starting with a certain index.We used the <a href=''http://www.pytables.org/moin''>Python bindings</a> to HDF. An added advantage of these bindings is that they will return the data back as numpy arrays (very fast numerical operations). </p>

<p>Here is the relevant code that deals with HDF only: <a href=''http://code.google.com/p/genetrack/source/browse/trunk/atlas/hdf.py''>hdf.py</a></p>

<p>The HDF schema is set up in a different module, but in the end it simply something like:</p>

<pre><code>class MySchema( IsDescription ):
    """
    Stores a triplet of float values for each index.
    """
    ix = IntCol  ( pos=1 )  # index
    wx = FloatCol( pos=2 )  # values on the W (forward) strand
    cx = FloatCol( pos=3 )  # value on the C (reverse) strand
    ax = FloatCol( pos=4 )  # weighted value on the combined W + C strands
</code></pre>
',69,0,'<p>In the <a rel="nofollow" href="http://code.google.com/p/genetrack/">GeneTrack</a> software we have used HDF to store values for each genomic base. Its main advantage over other storage systems was that it was able to return consecutive values with minimal overhead. </p>

<p>For example it is <em>extremely fast</em> (ms) in retrieving say 100,000 consecutive values starting with a certain index.We used the <a rel="nofollow" href="http://www.pytables.org/moin">Python bindings</a> to HDF. An added advantage of these bindings is that they will return the data back as numpy arrays (very fast numerical operations). </p>

<p>Here is the relevant code that deals with HDF only: <a rel="nofollow" href="http://code.google.com/p/genetrack/source/browse/trunk/atlas/hdf.py">hdf.py</a></p>

<p>The HDF schema is set up in a different module, but in the end it simply something like:</p>

<pre><code>class MySchema( IsDescription ):
    """
    Stores a triplet of float values for each index.
    """
    ix = IntCol  ( pos=1 )  # index
    wx = FloatCol( pos=2 )  # values on the W (forward) strand
    cx = FloatCol( pos=3 )  # value on the C (reverse) strand
    ax = FloatCol( pos=4 )  # weighted value on the combined W + C strands
</code></pre>
',1,1,0,69,2,0,'2010-02-26 20:01:39.187000',1,0,2);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-26 20:26:41.150000',0,'',74,0,0,'A: Using Hdf5 To Store  Bio-Data',0,4,'<p>What I do have is a netCDF-3 based Java application that I could show you.
NetCDF-3 is basically the same idea as HDF, but quite more limited as it cannot do compound datatypes among other limitations.</p>

<p>But here''s a small test code example to toy with:</p>

<p>package netCDF;</p>

<p>import java.io.File;
import ucar.ma2.<em>;
import ucar.nc2.</em>;
import java.io.IOException;
import java.util.ArrayList;</p>

<p>/**</p>

<hr />

<ul>
<li>@author Fernando Muñiz Fernandez</li>
<li>IBE, Institute of Evolutionary Biology (UPF-CSIC)</li>
<li>CEXS-UPF-PRBB</li>
</ul>

<hr />

<ul>
<li><p>THIS TO CREATE THE netCDF-3 GENOTYPE FILE
*/
public class CreateNetcdf {</p>

<p>public static NetcdfFileWriteable setDimsAndAttributes(Integer studyId, 
                                      String technology, 
                                      String description, 
                                      String strand, 
                                      int sampleSetSize,
                                      int markerSetSize) throws InvalidRangeException, IOException {</p>

<pre><code>///////////// CREATE netCDF-3 FILE ////////////
String genotypesFolder = "/media/data/genotypes";
File pathToStudy = new File(genotypesFolder+"/netCDF_test");
int gtSpan = constants.cNetCDF.Strides.STRIDE_GT;
int markerSpan = constants.cNetCDF.Strides.STRIDE_MARKER_NAME;
int sampleSpan = constants.cNetCDF.Strides.STRIDE_SAMPLE_NAME;

String matrixName = "prototype";
String writeFileName = pathToStudy+"/"+matrixName+".nc";
NetcdfFileWriteable ncfile = NetcdfFileWriteable.createNew(writeFileName, false);

// add dimensions
Dimension samplesDim = ncfile.addDimension("samples", sampleSetSize);
Dimension markersDim = ncfile.addDimension("markers", markerSetSize);
Dimension gtSpanDim = ncfile.addDimension("span", gtSpan);
ArrayList dims = new ArrayList();
dims.add(samplesDim);
dims.add(markersDim);
dims.add(gtSpanDim);

ArrayList markerGenotypeDims = new ArrayList();
markerGenotypeDims.add(markersDim);
markerGenotypeDims.add(markerSpan);

ArrayList markerPositionDim = new ArrayList();
markerPositionDim.add(markersDim);

ArrayList markerPropertyDim32 = new ArrayList();
markerPropertyDim32.add(markersDim);
markerPropertyDim32.add(32);

ArrayList markerPropertyDim16 = new ArrayList();
markerPropertyDim16.add(markersDim);
markerPropertyDim16.add(16);

ArrayList markerPropertyDim8 = new ArrayList();
markerPropertyDim8.add(markersDim);
markerPropertyDim8.add(8);

ArrayList markerPropertyDim2 = new ArrayList();
markerPropertyDim2.add(markersDim);
markerPropertyDim2.add(2);

ArrayList markerPropertyDim1 = new ArrayList();
markerPropertyDim1.add(markersDim);
markerPropertyDim1.add(1);

ArrayList sampleSetDims = new ArrayList();
sampleSetDims.add(samplesDim);
sampleSetDims.add(sampleSpan);

// Define Marker Variables
ncfile.addVariable("markerset", DataType.CHAR, markerGenotypeDims);
ncfile.addVariableAttribute("markerset", constants.cNetCDF.Attributes.LENGTH, markerSetSize);

ncfile.addVariable("marker_chromosome", DataType.CHAR, markerPropertyDim8);
ncfile.addVariable("marker_position", DataType.CHAR, markerPropertyDim32);
ncfile.addVariable("marker_position_int", DataType.INT, markerPositionDim);
ncfile.addVariable("marker_strand", DataType.CHAR, markerPropertyDim8);

ncfile.addVariable("marker_property_1", DataType.CHAR, markerPropertyDim1);
ncfile.addVariable("marker_property_2", DataType.CHAR, markerPropertyDim2);
ncfile.addVariable("marker_property_8", DataType.CHAR, markerPropertyDim8);
ncfile.addVariable("marker_property_16", DataType.CHAR, markerPropertyDim16);
ncfile.addVariable("marker_property_32", DataType.CHAR, markerPropertyDim32);

// Define Sample Variables
ncfile.addVariable("sampleset", DataType.CHAR, sampleSetDims);
ncfile.addVariableAttribute("sampleset", constants.cNetCDF.Attributes.LENGTH, sampleSetSize);

// Define Genotype Variables
ncfile.addVariable("genotypes", DataType.CHAR, dims);
ncfile.addVariableAttribute("genotypes", constants.cNetCDF.Attributes.GLOB_STRAND, "+/-");

// add global attributes
ncfile.addGlobalAttribute(constants.cNetCDF.Attributes.GLOB_STUDY, studyId);
ncfile.addGlobalAttribute(constants.cNetCDF.Attributes.GLOB_TECHNOLOGY, "INTERNAL");
ncfile.addGlobalAttribute(constants.cNetCDF.Attributes.GLOB_DESCRIPTION, "Matrix created by MOAPI through addition of 2 matrices");

return ncfile;
</code></pre>

<p>}
}</p></li>
</ul>

<p>Use the above in the following way:</p>

<p>package netCDF;</p>

<p>import java.util.List;
import ucar.ma2.<em>;
import ucar.nc2.</em>;
import java.io.IOException;</p>

<p>/**</p>

<hr />

<p>* @author Fernando Muñiz Fernandez
 * IBE, Institute of Evolutionary Biology (UPF-CSIC)
 * CEXS-UPF-PRBB</p>

<hr />

<p>* THIS TO GENERATE A netCDF-3 GENOTYPE DB
 */</p>

<p>public class TestWriteNetcdf {</p>

<pre><code>public static void main(String[] arg) throws InvalidRangeException, IOException {

    NetcdfFileWriteable ncfile = netCDF.CreateNetcdf.setDimsAndAttributes(0, 
                                                                      "INTERNAL", 
                                                                      "test in TestWriteNetcdf", 
                                                                      "+/-", 
                                                                      5,
                                                                      10);

    // create the file
    try {
        ncfile.create();
    } catch (IOException e) {
        System.err.println("ERROR creating file "+ncfile.getLocation()+"\n"+e);
    }


    ////////////// FILL''ER UP! ////////////////
    List&lt;Dimension&gt; dims = ncfile.getDimensions();
    Dimension samplesDim = dims.get(0);
    Dimension markersDim = dims.get(1);
    Dimension markerSpanDim = dims.get(2);

    ArrayChar charArray = new ArrayChar.D3(samplesDim.getLength(),markersDim.getLength(),markerSpanDim.getLength());
    int i,j;
    Index ima = charArray.getIndex();


    int method = 1;
    switch (method) {
        case 1: 
            // METHOD 1: Feed the complete genotype in one go
            for (i=0; i&lt;samplesDim.getLength(); i++) {
                for (j=0; j&lt;markersDim.getLength(); j++) {
                    char c = (char) ((char) j + 65);
                    String s = Character.toString(c) + Character.toString(c);
                    charArray.setString(ima.set(i,j,0),s);
                    System.out.println("SNP: "+i);
                }
            }
            break;
        case 2: 
            //METHOD 2: One snp at a time -&gt; feed in all samples
            for (i=0; i&lt;markersDim.getLength(); i++) {
                charArray.setString(ima.set(i,0), "s"+i+"I0");
                System.out.println("SNP: "+i);
            }
            break;
        case 3: 
            //METHOD 3: One sample at a time -&gt; feed in all snps
            break;
    }



    int[] offsetOrigin = new int[3]; //0,0
    try {
        ncfile.write("genotypes", offsetOrigin, charArray);
        //ncfile.write("genotype", origin, A);
    } catch (IOException e) {
        System.err.println("ERROR writing file");
    } catch (InvalidRangeException e) {
        e.printStackTrace();
    }

    // close the file
    try {
        ncfile.close();
    } catch (IOException e) {
        System.err.println("ERROR creating file "+ncfile.getLocation()+"\n"+e);
    }

}
</code></pre>

<p>}</p>
',69,0,'<p>What I do have is a netCDF-3 based Java application that I could show you.
NetCDF-3 is basically the same idea as HDF, but quite more limited as it cannot do compound datatypes among other limitations.</p>

<p>But here''s a small test code example to toy with:</p>

<p>package netCDF;</p>

<p>import java.io.File;
import ucar.ma2.<em>;
import ucar.nc2.</em>;
import java.io.IOException;
import java.util.ArrayList;</p>

<p>/**</p>

<hr>

<ul>
<li>@author Fernando Muñiz Fernandez</li>
<li>IBE, Institute of Evolutionary Biology (UPF-CSIC)</li>
<li>CEXS-UPF-PRBB</li>
</ul>

<hr>

<ul>
<li><p>THIS TO CREATE THE netCDF-3 GENOTYPE FILE
*/
public class CreateNetcdf {</p>

<p>public static NetcdfFileWriteable setDimsAndAttributes(Integer studyId, 
                                      String technology, 
                                      String description, 
                                      String strand, 
                                      int sampleSetSize,
                                      int markerSetSize) throws InvalidRangeException, IOException {</p>

<pre><code>///////////// CREATE netCDF-3 FILE ////////////
String genotypesFolder = "/media/data/genotypes";
File pathToStudy = new File(genotypesFolder+"/netCDF_test");
int gtSpan = constants.cNetCDF.Strides.STRIDE_GT;
int markerSpan = constants.cNetCDF.Strides.STRIDE_MARKER_NAME;
int sampleSpan = constants.cNetCDF.Strides.STRIDE_SAMPLE_NAME;

String matrixName = "prototype";
String writeFileName = pathToStudy+"/"+matrixName+".nc";
NetcdfFileWriteable ncfile = NetcdfFileWriteable.createNew(writeFileName, false);

// add dimensions
Dimension samplesDim = ncfile.addDimension("samples", sampleSetSize);
Dimension markersDim = ncfile.addDimension("markers", markerSetSize);
Dimension gtSpanDim = ncfile.addDimension("span", gtSpan);
ArrayList dims = new ArrayList();
dims.add(samplesDim);
dims.add(markersDim);
dims.add(gtSpanDim);

ArrayList markerGenotypeDims = new ArrayList();
markerGenotypeDims.add(markersDim);
markerGenotypeDims.add(markerSpan);

ArrayList markerPositionDim = new ArrayList();
markerPositionDim.add(markersDim);

ArrayList markerPropertyDim32 = new ArrayList();
markerPropertyDim32.add(markersDim);
markerPropertyDim32.add(32);

ArrayList markerPropertyDim16 = new ArrayList();
markerPropertyDim16.add(markersDim);
markerPropertyDim16.add(16);

ArrayList markerPropertyDim8 = new ArrayList();
markerPropertyDim8.add(markersDim);
markerPropertyDim8.add(8);

ArrayList markerPropertyDim2 = new ArrayList();
markerPropertyDim2.add(markersDim);
markerPropertyDim2.add(2);

ArrayList markerPropertyDim1 = new ArrayList();
markerPropertyDim1.add(markersDim);
markerPropertyDim1.add(1);

ArrayList sampleSetDims = new ArrayList();
sampleSetDims.add(samplesDim);
sampleSetDims.add(sampleSpan);

// Define Marker Variables
ncfile.addVariable("markerset", DataType.CHAR, markerGenotypeDims);
ncfile.addVariableAttribute("markerset", constants.cNetCDF.Attributes.LENGTH, markerSetSize);

ncfile.addVariable("marker_chromosome", DataType.CHAR, markerPropertyDim8);
ncfile.addVariable("marker_position", DataType.CHAR, markerPropertyDim32);
ncfile.addVariable("marker_position_int", DataType.INT, markerPositionDim);
ncfile.addVariable("marker_strand", DataType.CHAR, markerPropertyDim8);

ncfile.addVariable("marker_property_1", DataType.CHAR, markerPropertyDim1);
ncfile.addVariable("marker_property_2", DataType.CHAR, markerPropertyDim2);
ncfile.addVariable("marker_property_8", DataType.CHAR, markerPropertyDim8);
ncfile.addVariable("marker_property_16", DataType.CHAR, markerPropertyDim16);
ncfile.addVariable("marker_property_32", DataType.CHAR, markerPropertyDim32);

// Define Sample Variables
ncfile.addVariable("sampleset", DataType.CHAR, sampleSetDims);
ncfile.addVariableAttribute("sampleset", constants.cNetCDF.Attributes.LENGTH, sampleSetSize);

// Define Genotype Variables
ncfile.addVariable("genotypes", DataType.CHAR, dims);
ncfile.addVariableAttribute("genotypes", constants.cNetCDF.Attributes.GLOB_STRAND, "+/-");

// add global attributes
ncfile.addGlobalAttribute(constants.cNetCDF.Attributes.GLOB_STUDY, studyId);
ncfile.addGlobalAttribute(constants.cNetCDF.Attributes.GLOB_TECHNOLOGY, "INTERNAL");
ncfile.addGlobalAttribute(constants.cNetCDF.Attributes.GLOB_DESCRIPTION, "Matrix created by MOAPI through addition of 2 matrices");

return ncfile;
</code></pre>

<p>}
}</p></li>
</ul>

<p>Use the above in the following way:</p>

<p>package netCDF;</p>

<p>import java.util.List;
import ucar.ma2.<em>;
import ucar.nc2.</em>;
import java.io.IOException;</p>

<p>/**</p>

<hr>

<p>* @author Fernando Muñiz Fernandez
 * IBE, Institute of Evolutionary Biology (UPF-CSIC)
 * CEXS-UPF-PRBB</p>

<hr>

<p>* THIS TO GENERATE A netCDF-3 GENOTYPE DB
 */</p>

<p>public class TestWriteNetcdf {</p>

<pre><code>public static void main(String[] arg) throws InvalidRangeException, IOException {

    NetcdfFileWriteable ncfile = netCDF.CreateNetcdf.setDimsAndAttributes(0, 
                                                                      "INTERNAL", 
                                                                      "test in TestWriteNetcdf", 
                                                                      "+/-", 
                                                                      5,
                                                                      10);

    // create the file
    try {
        ncfile.create();
    } catch (IOException e) {
        System.err.println("ERROR creating file "+ncfile.getLocation()+"\n"+e);
    }


    ////////////// FILL''ER UP! ////////////////
    List&lt;Dimension&gt; dims = ncfile.getDimensions();
    Dimension samplesDim = dims.get(0);
    Dimension markersDim = dims.get(1);
    Dimension markerSpanDim = dims.get(2);

    ArrayChar charArray = new ArrayChar.D3(samplesDim.getLength(),markersDim.getLength(),markerSpanDim.getLength());
    int i,j;
    Index ima = charArray.getIndex();


    int method = 1;
    switch (method) {
        case 1: 
            // METHOD 1: Feed the complete genotype in one go
            for (i=0; i&lt;samplesDim.getLength(); i++) {
                for (j=0; j&lt;markersDim.getLength(); j++) {
                    char c = (char) ((char) j + 65);
                    String s = Character.toString(c) + Character.toString(c);
                    charArray.setString(ima.set(i,j,0),s);
                    System.out.println("SNP: "+i);
                }
            }
            break;
        case 2: 
            //METHOD 2: One snp at a time -&gt; feed in all samples
            for (i=0; i&lt;markersDim.getLength(); i++) {
                charArray.setString(ima.set(i,0), "s"+i+"I0");
                System.out.println("SNP: "+i);
            }
            break;
        case 3: 
            //METHOD 3: One sample at a time -&gt; feed in all snps
            break;
    }



    int[] offsetOrigin = new int[3]; //0,0
    try {
        ncfile.write("genotypes", offsetOrigin, charArray);
        //ncfile.write("genotype", origin, A);
    } catch (IOException e) {
        System.err.println("ERROR writing file");
    } catch (InvalidRangeException e) {
        e.printStackTrace();
    }

    // close the file
    try {
        ncfile.close();
    } catch (IOException e) {
        System.err.println("ERROR creating file "+ncfile.getLocation()+"\n"+e);
    }

}
</code></pre>

<p>}</p>
',1,1,0,69,42,0,'2010-02-26 20:26:41.150000',1,0,42);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-26 20:44:01.490000',0,'',75,0,0,'A: Site Use Guidelines',0,1,'<p>This is an excellent initiative: congratulations and thank you for setting it up!</p>

<p>I guess this site will be all the more useful as there are more contributers... So I guess that good questions for the administrator(s) of this site are: </p>

<ul>
<li>Do you have a plan for advertising this site/attracting new Users?</li>
<li>How can the Users help?</li>
</ul>
',1,0,'<p>This is an excellent initiative: congratulations and thank you for setting it up!</p>

<p>I guess this site will be all the more useful as there are more contributers... So I guess that good questions for the administrator(s) of this site are: </p>

<ul>
<li>Do you have a plan for advertising this site/attracting new Users?</li>
<li>How can the Users help?</li>
</ul>
',1,1,0,1,26,0,'2010-02-26 20:44:01.490000',1,0,26);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-26 21:37:25.810000',0,'compilation,taverna,plugin,maven,workflow',76,77,3,'Looking For A ''Hello World" Plugin For Taverna.',0,3,'<p>Hi all,
I''d like to create a very simple plugin for <a href=''http://www.taverna.org.uk/''>Taverna 2.0</a>, something very simple like like implementing a ''convertDnaToRna''. There is already some source code that can be found on the net e.g. Egon Willighagen''s code at <a href=''http://github.com/egonw/cdk-taverna''>http://github.com/egonw/cdk-taverna</a> but it requires to know <strong>Maven</strong> and.... I''m too <strong>lazy</strong> :-)</p>

<p>How can I implement this kind of simple plugin without maven ? ( I <em>just</em> want to compile, package &amp; create the right XML config files)</p>

<p>Thanks !</p>
',76,0,'<p>Hi all,
I''d like to create a very simple plugin for <a rel="nofollow" href="http://www.taverna.org.uk/">Taverna 2.0</a>, something very simple like like implementing a ''convertDnaToRna''. There is already some source code that can be found on the net e.g. Egon Willighagen''s code at <a rel="nofollow" href="http://github.com/egonw/cdk-taverna">http://github.com/egonw/cdk-taverna</a> but it requires to know <strong>Maven</strong> and.... I''m too <strong>lazy</strong> :-)</p>

<p>How can I implement this kind of simple plugin without maven ? ( I <em>just</em> want to compile, package &amp; create the right XML config files)</p>

<p>Thanks !</p>
',0,1,0,76,30,0,'2010-03-03 19:39:44.443000',1,7,30);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-26 21:49:18.627000',1,'bed,conversion,format',77,111,3,'How Do I Convert An Illumina Export File To Bed?',0,1,'<p>I have some illumina data generated from the latest version of the illumina pipeline (1.6.0) I need to convert my data into BED to view in ucsc genome browser.</p>

<p>This seems like it should be a fairly common task, however, I am unable to find any scripts to convert my data.</p>
',77,0,'<p>I have some illumina data generated from the latest version of the illumina pipeline (1.6.0) I need to convert my data into BED to view in ucsc genome browser.</p>

<p>This seems like it should be a fairly common task, however, I am unable to find any scripts to convert my data.</p>
',0,1,0,77,23,0,'2010-03-02 23:45:34.307000',1,14,44);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-26 22:15:18.997000',0,'',78,0,0,'A: How Do I Convert An Illumina Export File To Bed?',0,2,'<p>I found a script <a href=''http://mng.iop.kcl.ac.uk/site/node/378''>on another site</a>, Uses perl but I have not checked for correctness:</p>

<pre><code>#!/usr/bin/perl

use strict;
use warnings;
use IO::File;

my $filename = shift @ARGV;
die "Usage\n\tperl sorted2bed.pl s_X_sorted.txt &gt; s_X_sorted.bed\n" unless $filename;
chomp $filename;

my $fh = new IO::File;
$fh-&gt;open("&lt; $filename") or die "Can''t open file $filename for reading: $!";

my $count = 1;
while(my $line = &lt;$fh&gt;){
   warn "Line $count\n" if $count%1000 == 0;
   $count++;
   my @line = split "\t", $line;
   my $chr = $line[10];
   $chr =~ s/(.+)\.fa/$1/;
   #Illumina is 1-based, BED is 0-based
   my $start = $line[12]-1;
   my $read = $line[8];
   my $end = $start + length $read;
   my $strand = $line[13] eq ''F'' ? ''+'': ''-'';
   my $score = $line[15];
   my $bedline = "$chr\t$start\t$end\t$read\t$score\t$strand\n";
   print $bedline;
}
$fh-&gt;close;

warn "Done";
</code></pre>
',77,0,'<p>I found a script <a rel="nofollow" href="http://mng.iop.kcl.ac.uk/site/node/378">on another site</a>, Uses perl but I have not checked for correctness:</p>

<pre><code>#!/usr/bin/perl

use strict;
use warnings;
use IO::File;

my $filename = shift @ARGV;
die "Usage\n\tperl sorted2bed.pl s_X_sorted.txt &gt; s_X_sorted.bed\n" unless $filename;
chomp $filename;

my $fh = new IO::File;
$fh-&gt;open("&lt; $filename") or die "Can''t open file $filename for reading: $!";

my $count = 1;
while(my $line = &lt;$fh&gt;){
   warn "Line $count\n" if $count%1000 == 0;
   $count++;
   my @line = split "\t", $line;
   my $chr = $line[10];
   $chr =~ s/(.+)\.fa/$1/;
   #Illumina is 1-based, BED is 0-based
   my $start = $line[12]-1;
   my $read = $line[8];
   my $end = $start + length $read;
   my $strand = $line[13] eq ''F'' ? ''+'': ''-'';
   my $score = $line[15];
   my $bedline = "$chr\t$start\t$end\t$read\t$score\t$strand\n";
   print $bedline;
}
$fh-&gt;close;

warn "Done";
</code></pre>
',1,1,0,77,10,0,'2010-02-26 22:15:18.997000',1,0,10);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-26 22:49:35.150000',5,'general,make,pipeline,organization',79,276,23,'How To Organize A Pipeline Of Small Scripts Together?',0,4,'<p>In bioinformatics it is very common to end up with a lot of small scripts, each one with a different scope - plotting a chart, converting a file into another format, execute small operations - so it is very important to have a good way to clue them together, to define which should be executed before the others and so on.</p>

<p>How do you deal with the problem? Do you use makefiles, taverna workflows, batch scripts, or any other solution?</p>
',79,0,'<p>In bioinformatics it is very common to end up with a lot of small scripts, each one with a different scope - plotting a chart, converting a file into another format, execute small operations - so it is very important to have a good way to clue them together, to define which should be executed before the others and so on.</p>

<p>How do you deal with the problem? Do you use makefiles, taverna workflows, batch scripts, or any other solution?</p>
',0,1,2,79,23,0,'2010-03-02 23:45:01.483000',1,35,23);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-26 22:58:09.550000',0,'',80,0,0,'A: How To Organize A Pipeline Of Small Scripts Together?',0,4,'<p>I don''t have personal experience with this package but it is something that I plan to explore in the near future:</p>

<p><strong><a href=''http://code.google.com/p/ruffus/''>Ruffus </a></strong> a lightweight python module to run computational pipelines. </p>
',79,0,'<p>I don''t have personal experience with this package but it is something that I plan to explore in the near future:</p>

<p><strong><a rel="nofollow" href="http://code.google.com/p/ruffus/">Ruffus </a></strong> a lightweight python module to run computational pipelines. </p>
',1,1,0,79,2,0,'2010-02-26 22:58:09.550000',1,0,2);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-26 23:03:52.940000',0,'',81,0,0,'A: How To Organize A Pipeline Of Small Scripts Together?',0,5,'<p>My favorite way of defining pipelines is by writing Makefiles, about which you can find <a href=''http://software-carpentry.org/build.html''>a very good introduction</a> in Software Carpentry for Bioinformatics: <a href=''http://swc.scipy.org/lec/build.html''>http://swc.scipy.org/lec/build.html</a> .</p>

<p>Although they have been originally developed for compiling programs, Makefiles allow to define which operations are needed to create each file, with a declarative syntax that it is a bit old-style but still does its job. Each Makefile is composed of a set of rules, which define operations needed to calculate a file and that can be combined together to make a pipeline. Other advantages of makefiles are conditional execution of tasks, so you can stop the execution of a pipeline and get back to it later, without having to repeat calculations. However, one of the big disadvantages of Makefiles is its old syntax... in particular, rules are identified by the names of the files that they create, and there is no such thing as ''titles'' for rules, which make more tricky.</p>

<p>I think one of the best solutions would be to use <a href=''http://skam.sourceforge.net/skam-intro.html''>BioMake</a>, that allow to define tasks with titles that are not the name of the output files. To understand it better, look at <a href=''http://skam.sourceforge.net/skam-intro.html''>this example</a>: you see that each rule has a title and a series of parameters like its output, inputs, comments, etc.</p>

<p>Unfortunately, I can''t make biomake to run on my computer, as it requires very old dependencies and it is written in a very difficult perl. I have tried many alternatives and I think that <a href=''http://rake.rubyforge.org/files/doc/rational_rdoc.html''>rake</a> is the one that is more close to biomake, but unfortunately I don''t understand ruby''s syntax. </p>

<p>So, I am still looking for a good alternative... Maybe one day I will have to time to re-write BioMake in python :-)</p>
',79,0,'<p>My favorite way of defining pipelines is by writing Makefiles, about which you can find <a rel="nofollow" href="http://software-carpentry.org/build.html">a very good introduction</a> in Software Carpentry for Bioinformatics: <a rel="nofollow" href="http://swc.scipy.org/lec/build.html">http://swc.scipy.org/lec/build.html</a> .</p>

<p>Although they have been originally developed for compiling programs, Makefiles allow to define which operations are needed to create each file, with a declarative syntax that it is a bit old-style but still does its job. Each Makefile is composed of a set of rules, which define operations needed to calculate a file and that can be combined together to make a pipeline. Other advantages of makefiles are conditional execution of tasks, so you can stop the execution of a pipeline and get back to it later, without having to repeat calculations. However, one of the big disadvantages of Makefiles is its old syntax... in particular, rules are identified by the names of the files that they create, and there is no such thing as ''titles'' for rules, which make more tricky.</p>

<p>I think one of the best solutions would be to use <a rel="nofollow" href="http://skam.sourceforge.net/skam-intro.html">BioMake</a>, that allow to define tasks with titles that are not the name of the output files. To understand it better, look at <a rel="nofollow" href="http://skam.sourceforge.net/skam-intro.html">this example</a>: you see that each rule has a title and a series of parameters like its output, inputs, comments, etc.</p>

<p>Unfortunately, I can''t make biomake to run on my computer, as it requires very old dependencies and it is written in a very difficult perl. I have tried many alternatives and I think that <a rel="nofollow" href="http://rake.rubyforge.org/files/doc/rational_rdoc.html">rake</a> is the one that is more close to biomake, but unfortunately I don''t understand ruby''s syntax. </p>

<p>So, I am still looking for a good alternative... Maybe one day I will have to time to re-write BioMake in python :-)</p>
',1,1,0,79,23,0,'2010-02-26 23:03:52.940000',1,0,23);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-26 23:04:53.577000',0,'',82,0,0,'A: How To Organize A Pipeline Of Small Scripts Together?',0,3,'<p>Since I work a lot with Python, I usually write a wrapper method that embeds the external script/program, i.e. calls it, parses its output and returns the desired information. The ''glueing'' of several such methods then takes place within my Python code that calls all these wrappers. I guess that''s a very common thing to do.</p>

<p>Chris</p>
',79,0,'<p>Since I work a lot with Python, I usually write a wrapper method that embeds the external script/program, i.e. calls it, parses its output and returns the desired information. The ''glueing'' of several such methods then takes place within my Python code that calls all these wrappers. I guess that''s a very common thing to do.</p>

<p>Chris</p>
',1,1,0,79,47,0,'2010-02-26 23:04:53.577000',1,0,47);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-26 23:16:41.967000',0,'',83,0,0,'A: Where Can I Get The Secondary Structure Of A Protein?',0,4,'<p>If you have the PDB file then you can use the standard tool called DSSP , it is supposed to be the gold standard for obtaining secondary structure. In case you just have sequence then I personally prefer PSIPRED , it takes evolutionary information into account to predict the secondary structure . According to CASP evaluation it is one of the best secondary structure predictor available.</p>
',48,0,'<p>If you have the PDB file then you can use the standard tool called DSSP , it is supposed to be the gold standard for obtaining secondary structure. In case you just have sequence then I personally prefer PSIPRED , it takes evolutionary information into account to predict the secondary structure . According to CASP evaluation it is one of the best secondary structure predictor available.</p>
',1,1,0,48,7,0,'2010-02-26 23:16:41.967000',1,0,7);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-27 02:16:31.970000',0,'',84,0,0,'A: How To Organize A Pipeline Of Small Scripts Together?',0,3,'<p>My answer would be: don''t bother. I''ve often found that much of the scripts I write are never used again after the initial use. Therefore spending time using a complex framework that considers dependency between scripts is a waste because the results might be negative and you never visit the analysis again. Even if you do end up using the script multiple times a simple hacky bash script might be more than enough to meet the requirements.</p>

<p>There will however be the 1-2% of initial analyses that return a interesting result and therefore need to be expanded with more deeper investigation. I think this is the point to invest more time time in organising the project. For me I use Rake because it''s simple and allows me to write in the language I''m used to (Ruby).</p>

<p>Overall I think pragmatism is the important factor in computational biology. Just do enough to get the results you need and only invest more time when it''s necessary. There''s so many blind alleys in computational analysis of biological data it''s not worth investing too much of your time until it''s necessary.</p>
',79,0,'<p>My answer would be: don''t bother. I''ve often found that much of the scripts I write are never used again after the initial use. Therefore spending time using a complex framework that considers dependency between scripts is a waste because the results might be negative and you never visit the analysis again. Even if you do end up using the script multiple times a simple hacky bash script might be more than enough to meet the requirements.</p>

<p>There will however be the 1-2% of initial analyses that return a interesting result and therefore need to be expanded with more deeper investigation. I think this is the point to invest more time time in organising the project. For me I use Rake because it''s simple and allows me to write in the language I''m used to (Ruby).</p>

<p>Overall I think pragmatism is the important factor in computational biology. Just do enough to get the results you need and only invest more time when it''s necessary. There''s so many blind alleys in computational analysis of biological data it''s not worth investing too much of your time until it''s necessary.</p>
',1,1,0,79,53,0,'2010-02-27 02:16:31.970000',1,0,53);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-27 02:27:32.303000',0,'',85,0,0,'A: What Is The Best Way To Share Scripts Between Members Of A Lab?',0,2,'<p>My lab uses a network-attached storage unit which every Linux workstation mounts by NFS at startup. It was reasonably cheap -- a couple hundred dollars per TB. We also keep copies of public databases on there. We put data sets on there as we''re working on them, and also put the more important scripts in a Mercurial repositiory.</p>

<p>As Marcos and Istvan mentioned, a wiki integrated with your VCS would be wise, and Trac (trac.edgewall.org) is the obvious choice for that.</p>
',58,0,'<p>My lab uses a network-attached storage unit which every Linux workstation mounts by NFS at startup. It was reasonably cheap -- a couple hundred dollars per TB. We also keep copies of public databases on there. We put data sets on there as we''re working on them, and also put the more important scripts in a Mercurial repositiory.</p>

<p>As Marcos and Istvan mentioned, a wiki integrated with your VCS would be wise, and Trac trac.edgewall.org) is the obvious choice for that.</p>
',1,1,0,58,24,0,'2010-02-27 02:27:32.303000',1,0,24);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-27 02:44:04.547000',0,'',86,0,0,'A: What Is The Best Way To Share Scripts Between Members Of A Lab?',0,2,'<p>This might be useful .</p>

<p><a href=''http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000424''>A Quick Guide to Organizing Computational Biology Projects</a></p>
',58,0,'<p>This might be useful .</p>

<p><a rel="nofollow" href="http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000424">A Quick Guide to Organizing Computational Biology Projects</a></p>
',1,1,0,58,7,0,'2010-02-27 02:44:04.547000',1,0,7);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-02-28 03:00:56.140000',0,'',87,0,0,'A: Using Hdf5 To Store  Bio-Data',0,2,'<p>There is also a Perl binding to HDF5: PDL::IO::HDF5</p>

<p><a href=''http://search.cpan.org/~cerney/PDL-IO-HDF5-0.5/''>http://search.cpan.org/~cerney/PDL-IO-HDF5-0.5/</a>
This requires the Perl Data Language (PDL) package. The way, data-structures can be handled, sub-ranges of data can be defined  an data can be manipulated is actually very elegant in PDL such that computational code can profit from PDLs vectorized style of writing expressions.</p>

<p>The same is true for R and the hdf5 package: <a href=''http://cran.r-project.org/web/packages/hdf5/index.html''>http://cran.r-project.org/web/packages/hdf5/index.html</a></p>

<p>Code examples are in the package documentations of both, the R-hdf5 package documentation is quite little though.</p>

<p>Both of these language bindings might be a very efficient way to read and write HDF5 files.</p>

<p>There are also APIs in Fortran, Java, Python, Matlab, C, or C++. So it might make sense to select the language and define the type of data you wish to store first. </p>
',69,0,'<p>There is also a Perl binding to HDF5: PDL::IO::HDF5</p>

<p><a rel="nofollow" href="http://search.cpan.org/~cerney/PDL-IO-HDF5-0.5/">http://search.cpan.org/~cerney/PDL-IO-HDF5-0.5/</a>
This requires the Perl Data Language (PDL) package. The way, data-structures can be handled, sub-ranges of data can be defined  an data can be manipulated is actually very elegant in PDL such that computational code can profit from PDLs vectorized style of writing expressions.</p>

<p>The same is true for R and the hdf5 package: <a rel="nofollow" href="http://cran.r-project.org/web/packages/hdf5/index.html">http://cran.r-project.org/web/packages/hdf5/index.html</a></p>

<p>Code examples are in the package documentations of both, the R-hdf5 package documentation is quite little though.</p>

<p>Both of these language bindings might be a very efficient way to read and write HDF5 files.</p>

<p>There are also APIs in Fortran, Java, Python, Matlab, C, or C++. So it might make sense to select the language and define the type of data you wish to store first. </p>
',1,1,0,69,55,0,'2010-02-28 03:00:56.140000',1,0,55);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-03-01 19:51:12.300000',2,'general,agile,good,team',88,67,3,'Agile Programming For Bioinformaticians - Any Suggestions?',0,1,'<p>I am planning to prepare a talk for my workmates, to introduce them the basics of some agile programming methodology, which I think could give us good ideas to improve our working as a team.</p>

<p>My idea was to take inspiration from <a href=''http://www.extremeprogramming.org/rules/''>extreme programming</a> and explain the rules I like the most: use of <a href=''http://www.extremeprogramming.org/example/crcsim.html''>A7 cards to write tasks</a>, <a href=''http://www.extremeprogramming.org/rules/planninggame.html''>release planning</a> every 3 week, stand-up meeting every day, <a href=''http://www.extremeprogramming.org/rules/movepeople.html''>Move people around</a>, <a href=''http://www.extremeprogramming.org/rules/testfirst.html''>unit tests first</a>, <a href=''http://www.extremeprogramming.org/rules/pair.html''>pair programming</a> (at least introduce the concept), <a href=''http://www.extremeprogramming.org/rules/collective.html''>collective ownership</a>.</p>

<p>It is difficult for me to explain these rules as I don''t have much direct experience with, apart for few exceptions, and it is even more difficult because I will have to explain them to people who are not comfortable with programming and with software engineering in general.
However, I also think that I have to prepare this talk early and it will be much more difficult if I wait too much.</p>

<p>Do you have any experience with what I am talking about? Do you have any advice to give me, or can you recommend me a book or a practice that I could explain along with extreme programming?</p>
',88,0,'<p>I am planning to prepare a talk for my workmates, to introduce them the basics of some agile programming methodology, which I think could give us good ideas to improve our working as a team.</p>

<p>My idea was to take inspiration from <a rel="nofollow" href="http://www.extremeprogramming.org/rules/">extreme programming</a> and explain the rules I like the most: use of <a rel="nofollow" href="http://www.extremeprogramming.org/example/crcsim.html">A7 cards to write tasks</a>, <a rel="nofollow" href="http://www.extremeprogramming.org/rules/planninggame.html">release planning</a> every 3 week, stand-up meeting every day, <a rel="nofollow" href="http://www.extremeprogramming.org/rules/movepeople.html">Move people around</a>, <a rel="nofollow" href="http://www.extremeprogramming.org/rules/testfirst.html">unit tests first</a>, <a rel="nofollow" href="http://www.extremeprogramming.org/rules/pair.html">pair programming</a> (at least introduce the concept), <a rel="nofollow" href="http://www.extremeprogramming.org/rules/collective.html">collective ownership</a>.</p>

<p>It is difficult for me to explain these rules as I don''t have much direct experience with, apart for few exceptions, and it is even more difficult because I will have to explain them to people who are not comfortable with programming and with software engineering in general.
However, I also think that I have to prepare this talk early and it will be much more difficult if I wait too much.</p>

<p>Do you have any experience with what I am talking about? Do you have any advice to give me, or can you recommend me a book or a practice that I could explain along with extreme programming?</p>
',0,1,0,88,23,0,'2010-03-02 23:44:16.057000',1,21,23);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-03-01 20:06:55.973000',0,'',89,0,0,'A: Agile Programming For Bioinformaticians - Any Suggestions?',0,0,'<p>I think the approach is unsuited for individuals who are not comfortable with programming in general. There is a long way to go until someone becomes confident in their abilities. Before that this approach is not only ineffective, it might be even be detrimental.</p>

<p>Instead what helps most is transparency. Everyone needs to write code in a source code repository that can be viewed, commented and verified. People should become familiar with testing, code coverage, and continuous integration. </p>

<p>Something to read: <a href=''http://en.wikipedia.org/wiki/The_Mythical_Man-Month''>Mythical Man Month</a>.</p>
',88,0,'<p>I think the approach is unsuited for individuals who are not comfortable with programming in general. There is a long way to go until someone becomes confident in their abilities. Before that this approach is not only ineffective, it might be even be detrimental.</p>

<p>Instead what helps most is transparency. Everyone needs to write code in a source code repository that can be viewed, commented and verified. People should become familiar with testing, code coverage, and continuous integration. </p>

<p>Something to read: <a rel="nofollow" href="http://en.wikipedia.org/wiki/The_Mythical_Man-Month">Mythical Man Month</a>.</p>
',1,1,0,88,2,0,'2010-03-01 20:06:55.973000',1,0,2);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-03-01 20:26:09.937000',2,'sequence,biopython,python',90,112,5,'Computing The Reverse And Complement Of A Sequence With Biopython',0,1,'<p>An example that computes the reverse complement of a sequence with <a href=''http://biopython.org/wiki/Main_Page''>BioPython</a></p>

<pre><code>#
# Reverse complement example with BioPython
#

from Bio.Seq import Seq

# a separate function to reverse strings (or other iterables)
def rev(it):
    "Reverses an interable and returns it as a string"
    return ''''.join(reversed(it))

# create a Seq class instance
dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")

# original DNA
print type(dna)
print dna

# reverse complement DNA, returns a new sequence
print dna.reverse_complement()

# currently there is no direct way to just reverse a sequence
# we need to do a little extra work

rseq = rev(str(dna))
rdna = Seq(rseq)

# reversed sequence
print rdna

# to complement DNA, returns a new sequence
print dna.complement()
</code></pre>

<p>Produces the following output:</p>

<pre><code>&lt;class ''Bio.Seq.Seq''&gt;
ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT
GATAGCCCGTGGGAAAGTCGCCGGGTAATGTTACCGGTA
TACCGGTAACATTACCCGGCGACTTTCCCACGGGCTATC
</code></pre>
',90,0,'<p>An example that computes the reverse complement of a sequence with <a rel="nofollow" href="http://biopython.org/wiki/Main_Page">BioPython</a></p>

<pre><code>#
# Reverse complement example with BioPython
#

from Bio.Seq import Seq

# a separate function to reverse strings (or other iterables)
def rev(it):
    "Reverses an interable and returns it as a string"
    return ''''.join(reversed(it))

# create a Seq class instance
dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")

# original DNA
print type(dna)
print dna

# reverse complement DNA, returns a new sequence
print dna.reverse_complement()

# currently there is no direct way to just reverse a sequence
# we need to do a little extra work

rseq = rev(str(dna))
rdna = Seq(rseq)

# reversed sequence
print rdna

# to complement DNA, returns a new sequence
print dna.complement()
</code></pre>

<p>Produces the following output:</p>

<pre><code>&lt;class ''Bio.Seq.Seq''&gt;
ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT
GATAGCCCGTGGGAAAGTCGCCGGGTAATGTTACCGGTA
TACCGGTAACATTACCCGGCGACTTTCCCACGGGCTATC
</code></pre>
',0,1,0,90,10,0,'2010-03-01 22:11:52.887000',1,21,10);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-03-01 20:41:28.687000',0,'',91,0,0,'A: Agile Programming For Bioinformaticians - Any Suggestions?',0,2,'<p>I would suggest to have a look at <a href=''http://en.wikipedia.org/wiki/Scrum_(development)''>Scrum</a>, too. Certain parts would help not only bioinformations. For example estimating the time expenditure of tasks and the resulting burn down charts can be really helpful to see if something is stuck especially when working together on bigger projects.The daily scrum reports helps to meditate why who is doing what and offers a platform to discuss problems.</p>
',88,0,'<p>I would suggest to have a look at <a rel="nofollow" href="http://en.wikipedia.org/wiki/Scrum_(development)">Scrum</a>, too. Certain parts would help not only bioinformations. For example estimating the time expenditure of tasks and the resulting burn down charts can be really helpful to see if something is stuck especially when working together on bigger projects.The daily scrum reports helps to meditate why who is doing what and offers a platform to discuss problems.</p>
',1,1,0,88,39,0,'2010-03-01 20:41:28.687000',1,0,39);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-03-01 20:48:25.830000',0,'python,pygr,use,sequence',92,72,1,'Computing The Reverse And Complement Of A Sequence With Pygr',0,1,'<p>Computing the reverse complement with the <a href=''http://code.google.com/p/pygr/wiki/PygrDocumentation''>Pygr</a> bioinformatics framework:</p>

<pre><code>#
# Reverse complement example with pygr
#

from pygr.sequence import Sequence

# needs a separate function to reverse strings
def rev(it):
    "Reverses an interable and returns it as a string"
    return ''''.join(reversed(it))

# original sequence as as string
seq = ''ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG''

# create a Sequence class  instance named bobo
dna = Sequence(seq,''bobo'')

# sequence class'' type and content
print type(dna)
print dna

# the -operator reverse complements the DNA, returns a new sequence
print -dna

# to reverse the DNA, reverse the input data
rdna = Sequence( rev(seq),''bobo'')
print rdna

# to complement the DNA reverse complement, then reverse again
cseq = rev(str(-dna))
cdna = Sequence(cseq,''bobo'')

print cdna
</code></pre>

<p>Produces the output:</p>

<pre><code>&lt;class ''pygr.sequence.Sequence''&gt;
ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT
GATAGCCCGTGGGAAAGTCGCCGGGTAATGTTACCGGTA
TACCGGTAACATTACCCGGCGACTTTCCCACGGGCTATC
</code></pre>
',92,0,'<p>Computing the reverse complement with the <a rel="nofollow" href="http://code.google.com/p/pygr/wiki/PygrDocumentation">Pygr</a> bioinformatics framework:</p>

<pre><code>#
# Reverse complement example with pygr
#

from pygr.sequence import Sequence

# needs a separate function to reverse strings
def rev(it):
    "Reverses an interable and returns it as a string"
    return ''''.join(reversed(it))

# original sequence as as string
seq = ''ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG''

# create a Sequence class  instance named bobo
dna = Sequence(seq,''bobo'')

# sequence class'' type and content
print type(dna)
print dna

# the -operator reverse complements the DNA, returns a new sequence
print -dna

# to reverse the DNA, reverse the input data
rdna = Sequence( rev(seq),''bobo'')
print rdna

# to complement the DNA reverse complement, then reverse again
cseq = rev(str(-dna))
cdna = Sequence(cseq,''bobo'')

print cdna
</code></pre>

<p>Produces the output:</p>

<pre><code>&lt;class ''pygr.sequence.Sequence''&gt;
ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT
GATAGCCCGTGGGAAAGTCGCCGGGTAATGTTACCGGTA
TACCGGTAACATTACCCGGCGACTTTCCCACGGGCTATC
</code></pre>
',0,2,0,92,2,0,'2010-03-08 19:51:12.053000',1,7,10);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-03-01 20:51:52.190000',0,'',93,0,0,'A: Recommend Easy To Use Microarray Clustering Software',0,1,'<p>Possibly related:
<a href=''http://mmc.gnets.ncsu.edu/''>http://mmc.gnets.ncsu.edu/</a></p>
',5,0,'<p>Possibly related:
<a rel="nofollow" href="http://mmc.gnets.ncsu.edu/">http://mmc.gnets.ncsu.edu/</a></p>
',1,1,0,5,45,0,'2010-03-01 20:51:52.190000',1,0,45);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-03-01 22:07:15.283000',0,'',94,0,0,'A: Computing The Reverse And Complement Of A Sequence With Biopython',0,2,'<p>Wouldn''t it better to have a single question titled ''How to compute the reverse complement with python'' and put all the examples as different answers? Otherwise it seems a bit confusing..</p>
',90,0,'<p>Wouldn''t it better to have a single question titled ''How to compute the reverse complement with python'' and put all the examples as different answers? Otherwise it seems a bit confusing..</p>
',1,1,0,90,23,0,'2010-03-01 22:07:15.283000',1,0,23);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-03-01 22:32:26.107000',0,'',95,0,0,'A: How To Organize A Pipeline Of Small Scripts Together?',0,3,'<p>The most important thing for me has been keeping a README file at the top of each project directory, where I write down not just <em>how</em> to run the scripts, but <em>why</em> I wrote them in the first place -- coming back to a project after a several-month lull, it''s remarkable difficult to figure out what all the half-finished results mean without detailed notes.</p>

<p>That said:</p>

<ul>
<li><code>make</code> is pretty handy for simple pipelines that need to be re-run a lot</li>
<li>I''m also intrigued by <a href=''http://code.google.com/p/waf/''>waf</a> and <a href=''http://www.scons.org/''>scons</a>, since I use Python a lot</li>
<li>If a pipeline only takes a couple of minutes to run, and you only re-run it every few days, coercing it into a build system doesn''t really save time overall for that project</li>
<li>But once you''re used to working with a build system, the threshold where it pays off to use it on a new project drops dramatically</li>
</ul>
',79,0,'<p>The most important thing for me has been keeping a README file at the top of each project directory, where I write down not just <em>how</em> to run the scripts, but <em>why</em> I wrote them in the first place -- coming back to a project after a several-month lull, it''s remarkable difficult to figure out what all the half-finished results mean without detailed notes.</p>

<p>That said:</p>

<ul>
<li><code>make</code> is pretty handy for simple pipelines that need to be re-run a lot</li>
<li>I''m also intrigued by <a rel="nofollow" href="http://code.google.com/p/waf/">waf</a> and <a rel="nofollow" href="http://www.scons.org/">scons</a>, since I use Python a lot</li>
<li>If a pipeline only takes a couple of minutes to run, and you only re-run it every few days, coercing it into a build system doesn''t really save time overall for that project</li>
<li>But once you''re used to working with a build system, the threshold where it pays off to use it on a new project drops dramatically</li>
</ul>
',1,1,0,79,24,0,'2010-03-01 22:32:26.107000',1,0,24);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-03-01 22:52:22.510000',0,'',96,0,0,'A: Computing The Reverse And Complement Of A Sequence With Biopython',0,2,'<p>The Bio.Seq module provides two easy ways to get the complement and reverse complement from a sequence:</p>

<ul>
<li>If you have a string, use the functions <code>complement(dna)</code> and <code>reverse_complement(dna)</code></li>
<li>If you have a Seq object, use its methods with the same names: <code>dna.complement()</code> and <code>dna.reverse_complement</code></li>
</ul>

<p>To reverse a sequence, there is a function in the <code>Bio.SeqUtils</code> module called <code>reverse</code> which does what you would expect.</p>

<hr />

<p>(Sorry for going meta, but I don''t have commenting privileges yet. This can be deleted if the original post is edited.)</p>

<p>According to <a href=''http://meta.stackoverflow.com/questions/17845/etiquette-for-answering-your-own-question''>Meta Stack Overflow</a>, if you want to share the answer to a difficult question that''s poorly documented elsewhere online, you should post the question as a genuine one, and then submit your own answer separately. In theory, someone else may have an answer that''s better than yours, and this allows it to be voted to the top properly.</p>
',90,0,'<p>The Bio.Seq module provides two easy ways to get the complement and reverse complement from a sequence:</p>

<ul>
<li>If you have a string, use the functions <code>complement(dna)</code> and <code>reverse_complement(dna)</code></li>
<li>If you have a Seq object, use its methods with the same names: <code>dna.complement()</code> and <code>dna.reverse_complement</code></li>
</ul>

<p>To reverse a sequence, there is a function in the <code>Bio.SeqUtils</code> module called <code>reverse</code> which does what you would expect.</p>

<hr>

<p>(Sorry for going meta, but I don''t have commenting privileges yet. This can be deleted if the original post is edited.)</p>

<p>According to <a rel="nofollow" href="http://meta.stackoverflow.com/questions/17845/etiquette-for-answering-your-own-question">Meta Stack Overflow</a>, if you want to share the answer to a difficult question that''s poorly documented elsewhere online, you should post the question as a genuine one, and then submit your own answer separately. In theory, someone else may have an answer that''s better than yours, and this allows it to be voted to the top properly.</p>
',1,1,0,90,24,0,'2010-03-01 22:52:22.510000',1,0,24);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-03-01 23:05:13.920000',0,'',97,0,0,'A: Gene Id Conversion Tool',0,2,'<p><a href=''http://idconverter.bioinfo.cnio.es/''>http://idconverter.bioinfo.cnio.es/</a></p>

<p>Is another possible solution to this, although you might find this is not as up to date as you might like either.</p>
',22,0,'<p><a rel="nofollow" href="http://idconverter.bioinfo.cnio.es/">http://idconverter.bioinfo.cnio.es/</a></p>

<p>Is another possible solution to this, although you might find this is not as up to date as you might like either.</p>
',1,1,0,22,59,0,'2010-03-01 23:05:13.920000',1,0,59);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-03-02 01:51:25.980000',0,'',98,0,0,'A: Gene Id Conversion Tool',0,3,'<p>BioMart has already been mentioned. It can do much more than ID conversion but it is very useful for conversion purposes, it is regularly updated and you can select different genome builds and all kinds of genomic features. It seems to me that you wish to retrieve GeneIDs linked to Affymetrix IDs. To select these attributes in BioMart: go to the <a href=''http://www.biomart.org/biomart/martview''>Martview</a> page to start a new BioMart query.</p>

<p>Select attributes on the attribute page: The Ensembl GeneIDs and Transcript IDs are default. Ensembl GeneID and Affy IDs are under the "External" tab. Select your chip there.
To limit to those genes which are on the chip, use the Filters->Gene menue. You can limit the genes to those present on various platforms or your favourite set.</p>

<p>There is an URL button in biomart that allows to retrieve a URL for your query and to pass it on to others. Try this example:</p>

<p><a href=''http://www.biomart.org/biomart/martview?VIRTUALSCHEMANAME=default&amp;ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id|hsapiens_gene_ensembl.default.feature_page.embl|hsapiens_gene_ensembl.default.feature_page.affy_hg_u133a&amp;FILTERS=hsapiens_gene_ensembl.default.filters.with_affy_hg_u133a.only&amp;VISIBLEPANEL=resultspanel''>BioMart URL</a> URL, that should be a good starting point.</p>

<p>If you are interested in KEGG identifiers (Pathways, Genes), EC-numbers, etc. the  </p>

<p><a href=''http://www.genome.jp/kegg/kegg3.html''>KEGG Identifier page</a> could be handy, because the KEGG ids are not in BioMart as far as I know.</p>
',22,0,'<p>BioMart has already been mentioned. It can do much more than ID conversion but it is very useful for conversion purposes, it is regularly updated and you can select different genome builds and all kinds of genomic features. It seems to me that you wish to retrieve GeneIDs linked to Affymetrix IDs. To select these attributes in BioMart: go to the <a rel="nofollow" href="http://www.biomart.org/biomart/martview">Martview</a> page to start a new BioMart query.</p>

<p>Select attributes on the attribute page: The Ensembl GeneIDs and Transcript IDs are default. Ensembl GeneID and Affy IDs are under the "External" tab. Select your chip there.
To limit to those genes which are on the chip, use the Filters-&gt;Gene menue. You can limit the genes to those present on various platforms or your favourite set.</p>

<p>There is an URL button in biomart that allows to retrieve a URL for your query and to pass it on to others. Try this example:</p>

<p><a rel="nofollow" href="http://www.biomart.org/biomart/martview?VIRTUALSCHEMANAME=default&amp;ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id|hsapiens_gene_ensembl.default.feature_page.embl|hsapiens_gene_ensembl.default.feature_page.affy_hg_u133a&amp;FILTERS=hsapiens_gene_ensembl.default.filters.with_affy_hg_u133a.only&amp;VISIBLEPANEL=resultspanel">BioMart URL</a> URL, that should be a good starting point.</p>

<p>If you are interested in KEGG identifiers (Pathways, Genes), EC-numbers, etc. the  </p>

<p><a rel="nofollow" href="http://www.genome.jp/kegg/kegg3.html">KEGG Identifier page</a> could be handy, because the KEGG ids are not in BioMart as far as I know.</p>
',1,1,0,22,55,0,'2010-03-02 01:51:25.980000',1,0,55);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-03-02 02:01:20.047000',1,'interval,query,use,genomics',99,145,4,'Fast Interval Intersection Methodologies',0,2,'<p>Most genomic annotations are specified as intervals along the genome. </p>

<ul>
<li><a href=''http://books.google.com/books?id=NLngYyWFl_YC&amp;lpg=PA311&amp;ots=BwTtEE-jJ9&amp;dq=cormen%20interval%20tree&amp;pg=PA311#v=onepage&amp;q=&amp;f=false''>Interval trees</a> have been known to provide an efficient datastructure that allows for very fast overlap querying. </li>
<li><a href=''http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btl647v1''>Nested Containment Lists</a> have been proposed as an even faster alternative </li>
</ul>

<p>Provide code examples in your programming language that demonstrate the use of fast interval querying.</p>
',99,0,'<p>Most genomic annotations are specified as intervals along the genome. </p>

<ul>
<li><a rel="nofollow" href="http://books.google.com/books?id=NLngYyWFl_YC&amp;lpg=PA311&amp;ots=BwTtEE-jJ9&amp;dq=cormen%20interval%20tree&amp;pg=PA311#v=onepage&amp;q=&amp;f=false">Interval trees</a> have been known to provide an efficient datastructure that allows for very fast overlap querying. </li>
<li><a rel="nofollow" href="http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btl647v1">Nested Containment Lists</a> have been proposed as an even faster alternative </li>
</ul>

<p>Provide code examples in your programming language that demonstrate the use of fast interval querying.</p>
',0,1,1,99,23,0,'2010-03-02 23:42:32.173000',1,14,10);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2010-03-02 02:06:18.193000',0,'',100,0,0,'A: Fast Interval Intersection Methodologies',0,2,'<p>This code example generates 10,000 intervals then queries them for overlapping regions. <strong>Requires only the presence of Python.</strong></p>

<p>The code below requires the either the installation of the <a href=''http://bitbucket.org/james_taylor/bx-python/wiki/Home''>bx python</a> package or alternatively you may just download the <a href=''http://bitbucket.org/james_taylor/bx-python/raw/ebf9a4b352d3/lib/bx/intervals/operations/quicksect.py''>quicksect.py</a> module and place it next to the script itself:</p>

<pre><code>from random import randint, seed

# if you can install bx python then uncomment the line below
#
# from bx.intervals.operations.quicksect import IntervalNode

# otherwise just download the quickset module as shown above 
# and place it in next to your program
#
from quicksect import IntervalNode

# the span of the generated intervals
SPAN = 10

# the size of the genome
SIZE = 5*10**4

# the number of intervals
N = 10**4

def generate(x):
    "Generates random interval over a size and span"
    lo = randint(10000, SIZE)
    hi = lo + randint(1, SPAN)
    return (lo, hi)

def find(start, end, tree):
    "Returns a list with the overlapping intervals"
    out = []
    tree.intersect( start, end, lambda x: out.append(x) )
    return [ (x.start, x.end) for x in out ]

# use this to force both examples to generate the same data
seed(10)

# generate 10 thousand random intervals
data = map(generate, xrange(N))

# generate the intervals to query over
query = map(generate, xrange(10))

# start the root at the first element
start, end = data[0]
tree = IntervalNode( start, end )

# build an interval tree from the rest of the data
for start, end in data[1:]:
    tree = tree.insert( start, end )

for start, end in query:
    overlap = find(start, end , tree)
    print ''(%s, %s) -&gt; %s'' % (start, end, overlap)
</code></pre>

<p>Produces the output:</p>

<pre><code>(41901, 41903) -&gt; [(41894, 41902)]
(36981, 36987) -&gt; [(36981, 36984), (36973, 36982), (36978, 36987)]
(36338, 36339) -&gt; [(36337, 36347)]
(32741, 32748) -&gt; [(32738, 32742)]
(49864, 49872) -&gt; [(49859, 49865)]
(21475, 21477) -&gt; []
(29425, 29428) -&gt; [(29418, 29426), (29419, 29426)]
(29590, 29599) -&gt; [(29586, 29595), (29596, 29598)]
(12804, 12811) -&gt; [(12806, 12811), (12799, 12806), (12809, 12819)]
(30339, 30343) -&gt; [(30336, 30346), (30335, 30345), (30340, 30341)]
</code></pre>
',99,0,'<p>This code example generates 10,000 intervals then queries them for overlapping regions. <strong>Requires only the presence of Python.</strong></p>

<p>The code below requires the either the installation of the <a rel="nofollow" href="http://bitbucket.org/james_taylor/bx-python/wiki/Home">bx python</a> package or alternatively you may just download the <a rel="nofollow" href="http://bitbucket.org/james_taylor/bx-python/raw/ebf9a4b352d3/lib/bx/intervals/operations/quicksect.py">quicksect.py</a> module and place it next to the script itself:</p>

<pre><code>from random import randint, seed

# if you can install bx python then uncomment the line below
#
# from bx.intervals.operations.quicksect import IntervalNode

# otherwise just download the quickset module as shown above 
# and place it in next to your program
#
from quicksect import IntervalNode

# the span of the generated intervals
SPAN = 10

# the size of the genome
SIZE = 5*10**4

# the number of intervals
N = 10**4

def generate(x):
    "Generates random interval over a size and span"
    lo = randint(10000, SIZE)
    hi = lo + randint(1, SPAN)
    return (lo, hi)

def find(start, end, tree):
    "Returns a list with the overlapping intervals"
    out = []
    tree.intersect( start, end, lambda x: out.append(x) )
    return [ (x.start, x.end) for x in out ]

# use this to force both examples to generate the same data
seed(10)

# generate 10 thousand random intervals
data = map(generate, xrange(N))

# generate the intervals to query over
query = map(generate, xrange(10))

# start the root at the first element
start, end = data[0]
tree = IntervalNode( start, end )

# build an interval tree from the rest of the data
for start, end in data[1:]:
    tree = tree.insert( start, end )

for start, end in query:
    overlap = find(start, end , tree)
    print ''(%s, %s) -&gt; %s'' % (start, end, overlap)
</code></pre>

<p>Produces the output:</p>

<pre><code>(41901, 41903) -&gt; [(41894, 41902)]
(36981, 36987) -&gt; [(36981, 36984), (36973, 36982), (36978, 36987)]
(36338, 36339) -&gt; [(36337, 36347)]
(32741, 32748) -&gt; [(32738, 32742)]
(49864, 49872) -&gt; [(49859, 49865)]
(21475, 21477) -&gt; []
(29425, 29428) -&gt; [(29418, 29426), (29419, 29426)]
(29590, 29599) -&gt; [(29586, 29595), (29596, 29598)]
(12804, 12811) -&gt; [(12806, 12811), (12799, 12806), (12809, 12819)]
(30339, 30343) -&gt; [(30336, 30346), (30335, 30345), (30340, 30341)]
</code></pre>
',1,1,0,99,2,0,'2010-03-02 02:13:11.680000',1,0,2);
INSERT INTO "posts_post" VALUES(NULL,0.0,'2014-05-09 14:39:20.047000',0,'galaxy,motif',101,1,0,'How to find motifs with Galaxy?',0,0,'I have a few hundred yeast sequences (20-80bp long) and I want to find common motifs (conserved bases at certain indices) in them. I am using a Mac',101,0,'I have a few hundred yeast sequences (20-80bp long) and I want to find common motifs (conserved bases at certain indices) in them. I am using a Mac',0,1,0,101,1,0,'2014-05-09 14:39:20.047000',1,4,1);
CREATE TABLE "posts_replytoken" ("id" integer NOT NULL PRIMARY KEY, "date" datetime NOT NULL, "user_id" integer NOT NULL, "post_id" integer NOT NULL, "token" varchar(256) NOT NULL);
INSERT INTO "posts_replytoken" VALUES(1,'2014-04-29 15:30:13.832000',3,7,'58d43201e860670862b2cbb8d6177dc8');
INSERT INTO "posts_replytoken" VALUES(2,'2014-04-29 15:30:13.852000',3,8,'d76b18bacf7ab9403efa8a27b36b48d9');
INSERT INTO "posts_replytoken" VALUES(3,'2014-04-29 15:30:13.871000',3,9,'76d0cb5f7851acb2275dec01c17753a4');
INSERT INTO "posts_replytoken" VALUES(4,'2014-04-29 15:30:13.917000',10,11,'5a12c8be1838453c2fccce9028da087b');
INSERT INTO "posts_replytoken" VALUES(5,'2014-04-29 15:30:13.974000',12,14,'553be3c645911a8adc25388cda883e91');
INSERT INTO "posts_replytoken" VALUES(6,'2014-04-29 15:30:13.993000',12,15,'ec68a03d3eb79b24bfbe1efcf4171c76');
INSERT INTO "posts_replytoken" VALUES(7,'2014-04-29 15:30:14.013000',12,16,'90483e0fbf373bc6061823a7cac61ea6');
INSERT INTO "posts_replytoken" VALUES(8,'2014-04-29 15:30:14.032000',2,17,'d8056e9808ae7b32be0016b6204adf6d');
INSERT INTO "posts_replytoken" VALUES(9,'2014-04-29 15:30:14.051000',2,18,'efdab3ad7f7a289a745bd3d71a50015b');
INSERT INTO "posts_replytoken" VALUES(10,'2014-04-29 15:30:14.069000',2,19,'92b4f704fb00419c5b36787899130243');
INSERT INTO "posts_replytoken" VALUES(11,'2014-04-29 15:30:14.170000',16,23,'5d1a771b4c1912ca2895f20a6c2063b6');
INSERT INTO "posts_replytoken" VALUES(12,'2014-04-29 15:30:14.212000',14,25,'319df6b8e058bb0c1dd87be84a513654');
INSERT INTO "posts_replytoken" VALUES(13,'2014-04-29 15:30:14.234000',14,26,'896bcfb68e68d2f34dbea65ae1c71035');
INSERT INTO "posts_replytoken" VALUES(14,'2014-04-29 15:30:14.253000',16,27,'c47d61f4d27f0211b1c3022ee36b01b7');
INSERT INTO "posts_replytoken" VALUES(15,'2014-04-29 15:30:14.306000',20,29,'936baec406eed9c26e4e0caf0f2e6b1e');
INSERT INTO "posts_replytoken" VALUES(16,'2014-04-29 15:30:14.332000',2,30,'294deb383fb86bbd629379c3920d0725');
INSERT INTO "posts_replytoken" VALUES(17,'2014-04-29 15:30:14.375000',4,32,'ca90355c744d7e3d806cc778944b4235');
INSERT INTO "posts_replytoken" VALUES(18,'2014-04-29 15:30:14.437000',23,35,'8a6f4c1f812314c1e630b7d8f6be2811');
INSERT INTO "posts_replytoken" VALUES(19,'2014-04-29 15:30:14.457000',23,36,'019a948987c739155a6e6ff721ca8afc');
INSERT INTO "posts_replytoken" VALUES(20,'2014-04-29 15:30:14.482000',23,37,'1a68f57c15c32a4bdaa9d11a1a7dc90b');
INSERT INTO "posts_replytoken" VALUES(21,'2014-04-29 15:30:14.500000',23,38,'adf6fba6bf008c5f9be636e2ddd4ae5d');
INSERT INTO "posts_replytoken" VALUES(22,'2014-04-29 15:30:14.523000',16,39,'a5acd55d64d101085765a151dab13b2a');
INSERT INTO "posts_replytoken" VALUES(23,'2014-04-29 15:30:14.543000',3,40,'39e8bf6d9bfd261a53eca3d2a145df7f');
INSERT INTO "posts_replytoken" VALUES(24,'2014-04-29 15:30:14.586000',2,42,'bff8824e847f00173ff1bab7e8e5e976');
INSERT INTO "posts_replytoken" VALUES(25,'2014-04-29 15:30:14.629000',23,44,'7e533d86b1a8fd788e8dba32835847b4');
INSERT INTO "posts_replytoken" VALUES(26,'2014-04-29 15:30:14.648000',23,45,'38d2d47d4397b9b4cafe0e8ecedd7d72');
INSERT INTO "posts_replytoken" VALUES(27,'2014-04-29 15:30:14.689000',23,47,'0a6a60cc5e59d77a34f7a66643613ab2');
INSERT INTO "posts_replytoken" VALUES(28,'2014-04-29 15:30:14.729000',23,49,'7121feb4b1aa5c34a2b9112f10f252e2');
INSERT INTO "posts_replytoken" VALUES(29,'2014-04-29 15:30:14.751000',23,50,'3bab174f136d0aa1abd9e442e9487d55');
INSERT INTO "posts_replytoken" VALUES(30,'2014-04-29 15:30:14.791000',14,52,'13c887b3fcfce7aba7e6014ac5835ff0');
INSERT INTO "posts_replytoken" VALUES(31,'2014-04-29 15:30:14.831000',4,54,'d1e57b7541bdee21c2092009661ee8e5');
INSERT INTO "posts_replytoken" VALUES(32,'2014-04-29 15:30:14.850000',27,55,'d691fd0b75d5de98fad691ba3a27a653');
INSERT INTO "posts_replytoken" VALUES(33,'2014-04-29 15:30:14.892000',23,57,'f99b528c7799843808ff36f19d36621b');
INSERT INTO "posts_replytoken" VALUES(34,'2014-04-29 15:30:14.934000',23,59,'49a2df5601cb8ed59564f3d335d1bb34');
INSERT INTO "posts_replytoken" VALUES(35,'2014-04-29 15:30:14.953000',3,60,'b9730576b45d6ded5ee4bc2f1dc8ea54');
INSERT INTO "posts_replytoken" VALUES(36,'2014-04-29 15:30:14.972000',23,61,'500870618936ae72859fc2a1a6dc846d');
INSERT INTO "posts_replytoken" VALUES(37,'2014-04-29 15:30:14.992000',23,62,'a632e3202f7837d9d580e166a096176a');
INSERT INTO "posts_replytoken" VALUES(38,'2014-04-29 15:30:15.035000',23,63,'8693c0fa74972667167e60b9c300617e');
INSERT INTO "posts_replytoken" VALUES(39,'2014-04-29 15:30:15.054000',23,64,'02f793bbd89a552391f96f69114a888b');
INSERT INTO "posts_replytoken" VALUES(40,'2014-04-29 15:30:15.073000',23,65,'c43a17899c7d4d089003db532b3a8e2a');
INSERT INTO "posts_replytoken" VALUES(41,'2014-04-29 15:30:15.092000',2,66,'be3d2a90c8546168f90412ebf140eb02');
INSERT INTO "posts_replytoken" VALUES(42,'2014-04-29 15:30:15.112000',23,67,'264aa9e0a9797e926aacd8673a2c3756');
INSERT INTO "posts_replytoken" VALUES(43,'2014-04-29 15:30:15.134000',23,68,'6fceb91302ffb3465f81999a45d9c046');
INSERT INTO "posts_replytoken" VALUES(44,'2014-04-29 15:30:15.179000',23,70,'3f63e536fc85d8452d49ecfa3cec5824');
INSERT INTO "posts_replytoken" VALUES(45,'2014-04-29 15:30:15.201000',30,71,'1d9c15440b5f935f26ada22836c2fb29');
INSERT INTO "posts_replytoken" VALUES(46,'2014-04-29 15:30:15.224000',30,72,'cd5f5e3b77def753f1941042da0ac681');
INSERT INTO "posts_replytoken" VALUES(47,'2014-04-29 15:30:15.248000',30,73,'7712220c0a6a053b6d090bd1441dd012');
INSERT INTO "posts_replytoken" VALUES(48,'2014-04-29 15:30:15.284000',30,74,'80b2434e10ab57083331f38bea8ace21');
INSERT INTO "posts_replytoken" VALUES(49,'2014-04-29 15:30:15.306000',2,75,'cb32db62a6870b15e2681f1eb77b3808');
INSERT INTO "posts_replytoken" VALUES(50,'2014-04-29 15:30:15.372000',44,78,'3dae8388fbf78002be9ca4d180a88cec');
INSERT INTO "posts_replytoken" VALUES(51,'2014-04-29 15:30:15.411000',23,80,'89ff6f8e40efecce625afdbdd4f3bba0');
INSERT INTO "posts_replytoken" VALUES(52,'2014-04-29 15:30:15.454000',23,82,'e85c46580ef0c618d3345c55741fa9a7');
INSERT INTO "posts_replytoken" VALUES(53,'2014-04-29 15:30:15.473000',23,83,'876218bb0ff562fdb2d02735e0120495');
INSERT INTO "posts_replytoken" VALUES(54,'2014-04-29 15:30:15.493000',23,84,'6b65b8c1c5ea781ebc8b943d42fecaeb');
INSERT INTO "posts_replytoken" VALUES(55,'2014-04-29 15:30:15.514000',23,85,'56eeb0075d787a97358dcd3eded089ea');
INSERT INTO "posts_replytoken" VALUES(56,'2014-04-29 15:30:15.534000',23,86,'d278e262bed502bb2ec850f973c1b08d');
INSERT INTO "posts_replytoken" VALUES(57,'2014-04-29 15:30:15.557000',30,87,'e84f198025c3543d3d0646ce55a5511b');
INSERT INTO "posts_replytoken" VALUES(58,'2014-04-29 15:30:15.606000',23,89,'f0003fbf2765b0dd49f00ae353cbfab8');
INSERT INTO "posts_replytoken" VALUES(59,'2014-04-29 15:30:15.648000',23,91,'9b98ff409830dd448ed0011b1868c428');
INSERT INTO "posts_replytoken" VALUES(60,'2014-04-29 15:30:15.690000',2,93,'447e758d06f52403889c39ce4f66902c');
INSERT INTO "posts_replytoken" VALUES(61,'2014-04-29 15:30:15.709000',10,94,'57523b85c81c4c6bf2f97e03deac983f');
INSERT INTO "posts_replytoken" VALUES(62,'2014-04-29 15:30:15.734000',23,95,'05404dcb08c27d01b14d4279408e88ed');
INSERT INTO "posts_replytoken" VALUES(63,'2014-04-29 15:30:15.759000',10,96,'5e587d9f05dfc15226825a8de6abd311');
INSERT INTO "posts_replytoken" VALUES(64,'2014-04-29 15:30:15.781000',16,97,'8364111d4df17d5cc3a1935c2ac18cb9');
INSERT INTO "posts_replytoken" VALUES(65,'2014-04-29 15:30:15.805000',16,98,'5f1ae95e8a8d5253360ce9e2ebdf999c');
INSERT INTO "posts_replytoken" VALUES(66,'2014-04-29 15:30:15.856000',10,100,'b486a14648bb925dc6e43e52ea2e9108');
INSERT INTO "posts_replytoken" VALUES(67,'2014-05-09 14:39:13.584000',16,98,'8f0d04ce');
INSERT INTO "posts_replytoken" VALUES(68,'2014-05-09 14:39:13.639000',10,100,'1fd599a1');
INSERT INTO "posts_replytoken" VALUES(69,'2014-05-09 14:39:20.151000',2,101,'9a7e846c');
CREATE TABLE "posts_emailentry" ("id" integer NOT NULL PRIMARY KEY, "post_id" integer NULL, "text" text NOT NULL, "creation_date" datetime NOT NULL, "sent_at" datetime NULL, "status" integer NOT NULL);
CREATE TABLE "posts_emailsub" ("id" integer NOT NULL PRIMARY KEY, "email" varchar(75) NOT NULL, "status" integer NOT NULL);
CREATE TABLE "badges_badge" ("count" integer NOT NULL, "name" varchar(50) NOT NULL, "icon" varchar(250) NOT NULL, "unique" bool NOT NULL, "type" integer NOT NULL, "id" integer PRIMARY KEY, "desc" varchar(200) NOT NULL);
INSERT INTO "badges_badge" VALUES(0,'Autobiographer','fa fa-bullhorn',0,0,1,'has more than 80 characters in the information field of the user''s profile');
INSERT INTO "badges_badge" VALUES(0,'Student','fa fa-certificate',0,0,2,'asked a question with at least 3 up-votes');
INSERT INTO "badges_badge" VALUES(0,'Teacher','fa fa-smile-o',0,0,3,'created an answer with at least 3 up-votes');
INSERT INTO "badges_badge" VALUES(0,'Commentator','fa fa-comment',0,0,4,'created a comment with at least 3 up-votes');
INSERT INTO "badges_badge" VALUES(0,'Supporter','fa fa-thumbs-up',0,1,5,'voted at least 25 times');
INSERT INTO "badges_badge" VALUES(0,'Scholar','fa fa-check-circle-o',0,0,6,'created an answer that has been accepted');
INSERT INTO "badges_badge" VALUES(0,'Voter','fa fa-thumbs-o-up',0,0,7,'voted more than 100 times');
INSERT INTO "badges_badge" VALUES(0,'Centurion','fa fa-bolt',0,1,8,'created 100 posts');
INSERT INTO "badges_badge" VALUES(0,'Cylon','fa fa-rocket',0,2,9,'received 1,000 up votes');
INSERT INTO "badges_badge" VALUES(0,'Rising Star','fa fa-star',0,2,10,'created 50 posts within first three months of joining');
INSERT INTO "badges_badge" VALUES(0,'Guru','fa fa-beer',0,1,11,'received more than 100 upvotes');
INSERT INTO "badges_badge" VALUES(0,'Popular Question','fa fa-eye',0,2,12,'created a question with more than 1,000 views');
INSERT INTO "badges_badge" VALUES(0,'Epic Question','fa fa-bullseye',0,2,13,'created a question with more than 10,000 views');
INSERT INTO "badges_badge" VALUES(0,'Oracle','fa fa-sun-o',0,2,14,'created more than 1,000 posts (questions + answers + comments)');
INSERT INTO "badges_badge" VALUES(0,'Pundit','fa fa-comments-o',0,1,15,'created a comment with more than 10 votes');
INSERT INTO "badges_badge" VALUES(0,'Good Answer','fa fa-pencil-square-o',0,0,16,'created an answer that was upvoted at least 5 times');
INSERT INTO "badges_badge" VALUES(0,'Good Question','fa fa-question',0,0,17,'asked a question that was upvoted at least 5 times');
INSERT INTO "badges_badge" VALUES(0,'Prophet','fa fa-pagelines',0,0,18,'created a post with more than 20 followers');
INSERT INTO "badges_badge" VALUES(0,'Librarian','fa fa-bookmark-o',0,0,19,'created a post with more than 10 bookmarks');
INSERT INTO "badges_badge" VALUES(0,'Great Question','fa fa-fire',0,1,20,'created a question with more than 5,000 views');
INSERT INTO "badges_badge" VALUES(0,'Gold Standard','fa fa-bookmark',0,2,21,'created a post with more than 25 bookmarks');
INSERT INTO "badges_badge" VALUES(0,'Appreciated','fa fa-heart',0,1,22,'created a post with more than 5 votes');
CREATE TABLE "badges_award" ("date" datetime NOT NULL, "badge_id" integer NOT NULL, "user_id" integer NOT NULL, "id" integer PRIMARY KEY, "context" varchar(1000) NOT NULL);
CREATE TABLE "planet_blogpost" ("id" integer NOT NULL PRIMARY KEY, "blog_id" integer NOT NULL, "uid" varchar(200) NOT NULL, "title" varchar(200) NOT NULL, "content" text NOT NULL, "html" text NOT NULL, "creation_date" datetime NOT NULL, "insert_date" datetime NULL, "published" bool NOT NULL, "link" varchar(200) NOT NULL);
INSERT INTO "planet_blogpost" VALUES(1,1,'tag:blogger.com,1999:blog-8959227089815463704.post-6750884431982311373','Geographic population structure analysis of worldwide human populations infers their biogeographical origins : Nature Communications : Nature Publishing Group','  There''s no stopping a cheesy geneticist when he puts his mind to making a research idea stick ... like coming up with acronyms like GPS to determine biogeographical origin ...   http://www.nature.com/ncomms/2014/140429/ncomms4513/full/ncomms4513.html  ','  There''s no stopping a cheesy geneticist when he puts his mind to making a research idea stick ... like coming up with acronyms like GPS to determine biogeographical origin ...   http://www.nature.com/ncomms/2014/140429/ncomms4513/full/ncomms4513.html  ','2014-05-02 00:00:00','2014-05-09 14:47:06.262000',0,'http://feedproxy.google.com/~r/MyWeblogOnBioinformaticsGenomeScienceNextGenerationSequencing/~3/ALuC5xFhtc8/geographic-population-structure.html');
INSERT INTO "planet_blogpost" VALUES(2,1,'tag:blogger.com,1999:blog-8959227089815463704.post-894146776514373738','Gut check: Microbes in our stomachs may be making us miserable - Salon.com','Seems like gut microbiome is getting more and more air time in news.   http://www.salon.com/2014/04/28/gut_check_microbes_in_your_stomach_may_be_making_you_miserable_partner/?utm_source=facebook&amp;utm_medium=socialflow    ','Seems like gut microbiome is getting more and more air time in news.   http://www.salon.com/2014/04/28/gut_check_microbes_in_your_stomach_may_be_making_you_miserable_partner/?utm_source=facebook&amp;utm_medium=socialflow    ','2014-04-28 00:00:00','2014-05-09 14:47:06.269000',0,'http://feedproxy.google.com/~r/MyWeblogOnBioinformaticsGenomeScienceNextGenerationSequencing/~3/h_9na9QwXYs/gut-check-microbes-in-our-stomachs-may.html');
INSERT INTO "planet_blogpost" VALUES(3,1,'tag:blogger.com,1999:blog-8959227089815463704.post-2952401788432414585','Fwd: Welcome to the Google Genomics Preview','Sweet!!---------- Forwarded message ----------  Welcome to the Google Genomics Preview! You''ve been approved for early access to the API.       The goal of the Genomics API is to encourage interoperability and build a foundation to store, process, search, analyze and share tens of petabytes of genomic data.       We''ve loaded sample data from public BAM files:      * The complete 1000 Genomes Project      * Selections from the Personal Genome Project      How to get started:      * Follow the instructions in the developer documentation       * Try the sample genome browser which calls the API      * Try out the other open source examples -- an R script, Python MapReduce, and a Java file-based implementation      * Write your own code to call the API and explore new uses      This is only the beginning. Your feedback will be essential to make the API useful. Please submit feature requests, bugs and suggestions on our GitHub page.       Thank you for being part of the first wave. If you''d rather join with a different email address (Gmail or Google Apps domain), please fill out the request form with that address too, and we''ll grant access soon. Thank you for your interest!      Sincerely,      The Google Genomics team      ','Sweet!!---------- Forwarded message ----------  Welcome to the Google Genomics Preview! You''ve been approved for early access to the API.       The goal of the Genomics API is to encourage interoperability and build a foundation to store, process, search, analyze and share tens of petabytes of genomic data.       We''ve loaded sample data from public BAM files:      * The complete 1000 Genomes Project      * Selections from the Personal Genome Project      How to get started:      * Follow the instructions in the developer documentation       * Try the sample genome browser which calls the API      * Try out the other open source examples -- an R script, Python MapReduce, and a Java file-based implementation      * Write your own code to call the API and explore new uses      This is only the beginning. Your feedback will be essential to make the API useful. Please submit feature requests, bugs and suggestions on our GitHub page.       Thank you for being part of the first wave. If you''d rather join with a different email address (Gmail or Google Apps domain), please fill out the request form with that address too, and we''ll grant access soon. Thank you for your interest!      Sincerely,      The Google Genomics team      ','2014-04-20 00:00:00','2014-05-09 14:47:06.281000',0,'http://feedproxy.google.com/~r/MyWeblogOnBioinformaticsGenomeScienceNextGenerationSequencing/~3/cujIVDCBhYo/fwd-welcome-to-google-genomics-preview.html');
INSERT INTO "planet_blogpost" VALUES(4,2,'tag:blogger.com,1999:blog-36768584.post-5086008806446578077','NGS Saves A Young Life','One of the most electrifying talks at AGBT this year was given by Joe DeRisi of UCSF, who gave a brief intro on the difficulty of diagnosing the root cause of encephalitis (as it can be autoimmune, viral, protozoal, bacterial and probably a few other causes) and then ran down a gripping case history which seemed straight out of House.Read more »','One of the most electrifying talks at AGBT this year was given by Joe DeRisi of UCSF, who gave a brief intro on the difficulty of diagnosing the root cause of encephalitis (as it can be autoimmune, viral, protozoal, bacterial and probably a few other causes) and then ran down a gripping case history which seemed straight out of House.Read more »','2014-02-26 00:00:00','2014-05-09 14:47:07.033000',0,'http://omicsomics.blogspot.com/2014/02/ngs-saves-young-life.html');
INSERT INTO "planet_blogpost" VALUES(5,2,'tag:blogger.com,1999:blog-36768584.post-2256957260103276078','A Sunset for Draft Genomes?','The sun set during AGBT 2014 for a final time over a week ago.  The posters have long been down, and perhaps the liver enzyme levels of the attendees are now down to normal as well.  This year’s conference underscored a possibility that was suggested last year: that the era of the poorly connected, low quality draft genome is headed for the sunset as wellRead more »','The sun set during AGBT 2014 for a final time over a week ago.  The posters have long been down, and perhaps the liver enzyme levels of the attendees are now down to normal as well.  This year’s conference underscored a possibility that was suggested last year: that the era of the poorly connected, low quality draft genome is headed for the sunset as wellRead more »','2014-02-25 00:00:00','2014-05-09 14:47:07.051000',0,'http://omicsomics.blogspot.com/2014/02/a-sunset-for-draft-genomes.html');
INSERT INTO "planet_blogpost" VALUES(6,2,'tag:blogger.com,1999:blog-36768584.post-6988867569351222907','How will you deal with GRCh38?','I was foolishly attempting to catch up with Twitter last night during Valerie Schneider''s AGBT talk last night on the new human reference, GRCh38. After all, my personal answer to my title is nothing, because this isn''t a field I work in.  But Dr. Schneider is a very good speaker and I could not help but have my attention pulled in.  While clearly not the final word on a human reference, this new edition fixes many gaps, expands the coverage of highly polymorphic regions, and even models the difficult to assemble centromeres.  Better assembly, combined with emerging tools to better handle those complex regions via graph representations, means better mapping send better variant calls.So, a significant advance, but a bit unpleasant one if you are in the space.  You now have several ugly options before you with regard to your prior data mapped to an earlier reference.The do nothing option must appeal to some. Forgo the advantages of the new reference and just stick to the old. Perhaps start new projects on the new one, leading to a cacophony of internal tools dealing with different versions, with an ongoing risk of mismatched results. Also, cross your fingers that none of changes might be revised if analyzed against the new reference.  Perhaps this route will be rationalized as healthy procrastination until a well-vetted set of graph-aware mappers exist, but once you start putting-off it is hard to stop doing so. The other pole would be to embrace the new reference whole-heartedly and realign all the old data against the new reference. After burning a lot of compute cycles and storage space running in place, spend a lot of time reconciling old and new results. Then decide whether to ditch all your old alignments, or suffer an even larger storage burden.A tempting shortcut would be to just remap alignments and variants by the known relationships between the two references. After all, the vast majority of the results will simply shift coordinates a bit, but with no other effects.  In theory, one could estimate all the map regions that are now suspect and simply realign the reads which map to those regions, plus attempt to place reads that previously failed to map. Again reconciliation of results, but on a much reduced scale.None would seem particularly appealing options. Perhaps that latter route will be a growth industry of new tools acting on BAM, CRAM or VCF which themselves will provide a morass of competing claims of accuracy, efficiency and speed. Doesn''t make me at all in a hurry to leave a cozy world of haploid genomes that are often finished by a simple pipeline!','I was foolishly attempting to catch up with Twitter last night during Valerie Schneider''s AGBT talk last night on the new human reference, GRCh38. After all, my personal answer to my title is nothing, because this isn''t a field I work in.  But Dr. Schneider is a very good speaker and I could not help but have my attention pulled in.  While clearly not the final word on a human reference, this new edition fixes many gaps, expands the coverage of highly polymorphic regions, and even models the difficult to assemble centromeres.  Better assembly, combined with emerging tools to better handle those complex regions via graph representations, means better mapping send better variant calls.So, a significant advance, but a bit unpleasant one if you are in the space.  You now have several ugly options before you with regard to your prior data mapped to an earlier reference.The do nothing option must appeal to some. Forgo the advantages of the new reference and just stick to the old. Perhaps start new projects on the new one, leading to a cacophony of internal tools dealing with different versions, with an ongoing risk of mismatched results. Also, cross your fingers that none of changes might be revised if analyzed against the new reference.  Perhaps this route will be rationalized as healthy procrastination until a well-vetted set of graph-aware mappers exist, but once you start putting-off it is hard to stop doing so. The other pole would be to embrace the new reference whole-heartedly and realign all the old data against the new reference. After burning a lot of compute cycles and storage space running in place, spend a lot of time reconciling old and new results. Then decide whether to ditch all your old alignments, or suffer an even larger storage burden.A tempting shortcut would be to just remap alignments and variants by the known relationships between the two references. After all, the vast majority of the results will simply shift coordinates a bit, but with no other effects.  In theory, one could estimate all the map regions that are now suspect and simply realign the reads which map to those regions, plus attempt to place reads that previously failed to map. Again reconciliation of results, but on a much reduced scale.None would seem particularly appealing options. Perhaps that latter route will be a growth industry of new tools acting on BAM, CRAM or VCF which themselves will provide a morass of competing claims of accuracy, efficiency and speed. Doesn''t make me at all in a hurry to leave a cozy world of haploid genomes that are often finished by a simple pipeline!','2014-02-13 00:00:00','2014-05-09 14:47:07.059000',0,'http://omicsomics.blogspot.com/2014/02/how-will-you-deal-with-grch38.html');
INSERT INTO "planet_blogpost" VALUES(7,3,'tag:blogger.com,1999:blog-14688252.post-1413366006721621665','Parallelizing #RStats using #make','In the current post, I''ll show how to use R as the main SHELL of GNU-Make instead of using a classical linux shell like ''bash''. Why would you do this ? awesomeness
Make-based workflow management
Make-based execution with --jobs.  GNU make knows how to execute several recipes at once. Normally, make will execute only one recipe at a time, waiting for it to finish before executing the next. However','In the current post, I''ll show how to use R as the main SHELL of GNU-Make instead of using a classical linux shell like ''bash''. Why would you do this ? awesomeness
Make-based workflow management
Make-based execution with --jobs.  GNU make knows how to execute several recipes at once. Normally, make will execute only one recipe at a time, waiting for it to finish before executing the next. However','2014-01-30 00:00:00','2014-05-09 14:47:07.839000',0,'http://plindenbaum.blogspot.com/2014/01/parallelizing-rstats-using-make.html');
INSERT INTO "planet_blogpost" VALUES(8,3,'tag:blogger.com,1999:blog-14688252.post-4555212046957957643','Mapping the UCSC/Web-Sequences to a world map.','People at the UCSC have recently released a new track for the GenomeBrowser
We BLATted the Internet! The DNA sequences from 40 billion webpages mapped to hg19 and other species: http://t.co/5XAsFCguE2/ UCSC Genome Browser (@GenomeBrowser) January 23, 2014"We''re pleased to announce the release of the Web Sequences track on the UCSC Genome Browser. This track, produced in collaboration with','People at the UCSC have recently released a new track for the GenomeBrowser
We BLATted the Internet! The DNA sequences from 40 billion webpages mapped to hg19 and other species: http://t.co/5XAsFCguE2/ UCSC Genome Browser (@GenomeBrowser) January 23, 2014"We''re pleased to announce the release of the Web Sequences track on the UCSC Genome Browser. This track, produced in collaboration with','2014-01-29 00:00:00','2014-05-09 14:47:07.846000',0,'http://plindenbaum.blogspot.com/2014/01/mapping-ucscweb-sequences-to-world-map.html');
INSERT INTO "planet_blogpost" VALUES(9,3,'tag:blogger.com,1999:blog-14688252.post-842379474834141569','(video) #renabigo2013 "Make &amp; Bioinformatics:everything but #usegalaxy"','(In French): "Mon Make à moi": parallélisations, workflows et pipelines pour le NGS, tout sauf Galaxy [34:50 mn]
Les workflows en Bio-Informatique
Centre Inria de Rennes - IRISA
29 Novembre 2013','(In French): "Mon Make à moi": parallélisations, workflows et pipelines pour le NGS, tout sauf Galaxy [34:50 mn]
Les workflows en Bio-Informatique
Centre Inria de Rennes - IRISA
29 Novembre 2013','2013-12-20 00:00:00','2014-05-09 14:47:07.851000',0,'http://plindenbaum.blogspot.com/2013/12/video-renabigo2013-make.html');
INSERT INTO "planet_blogpost" VALUES(10,4,'http://nsaunders.wordpress.com/?p=3727','When is db=all not db=all? When you use Entrez ELink.','Just a brief technical note. I figured that for a given compound in PubChem, it would be interesting to know whether that compound had been used in a high-throughput experiment, which you might find in GEO. Very easy using the E-utilities, as implemented in the R package rentrez: Browsing the rentrez documentation, I note that […]','Just a brief technical note. I figured that for a given compound in PubChem, it would be interesting to know whether that compound had been used in a high-throughput experiment, which you might find in GEO. Very easy using the E-utilities, as implemented in the R package rentrez: Browsing the rentrez documentation, I note that […]','2014-04-29 00:00:00','2014-05-09 14:47:08.508000',0,'http://nsaunders.wordpress.com/2014/04/30/when-is-dball-not-dball-when-you-use-entrez-elink/');
INSERT INTO "planet_blogpost" VALUES(11,4,'http://nsaunders.wordpress.com/?p=3724','On the road: CSS and eResearch Conference 2014','Next week I’ll be in Melbourne for one of my favourite meetings, the annual Computational and Simulation Sciences and eResearch Conference. The main reason for my visit is the Bioinformatics FOAM workshop. Day 1 (March 27) is not advertised since it is an internal CSIRO day, but I’ll be presenting a talk titled “SQL, noSQL […]','Next week I’ll be in Melbourne for one of my favourite meetings, the annual Computational and Simulation Sciences and eResearch Conference. The main reason for my visit is the Bioinformatics FOAM workshop. Day 1 (March 27) is not advertised since it is an internal CSIRO day, but I’ll be presenting a talk titled “SQL, noSQL […]','2014-03-20 00:00:00','2014-05-09 14:47:08.515000',0,'http://nsaunders.wordpress.com/2014/03/20/on-the-road-css-and-eresearch-conference-2014/');
INSERT INTO "planet_blogpost" VALUES(12,4,'http://nsaunders.wordpress.com/?p=3675','“Advance” access and DOIs: what’s the problem?','When I arrive at work, the first task for the day is “check feeds”. If I’m lucky, in the “journal TOCs” category, there will be an abstract that looks interesting, like this one on the left (click for larger version). Sometimes, the title is a direct link to the article at the journal website. Often […]','When I arrive at work, the first task for the day is “check feeds”. If I’m lucky, in the “journal TOCs” category, there will be an abstract that looks interesting, like this one on the left (click for larger version). Sometimes, the title is a direct link to the article at the journal website. Often […]','2014-03-09 00:00:00','2014-05-09 14:47:08.521000',0,'http://nsaunders.wordpress.com/2014/03/10/advance-access-and-dois-whats-the-problem/');
INSERT INTO "planet_blogpost" VALUES(13,5,'http://bcbio.wordpress.com/?p=570','Improving reproducibility and installation of genomic analysis pipelines with Docker','Motivation bcbio-nextgen is a community developed, best-practice pipeline for genomic data processing, performing automated variant calling and RNA-seq analyses from high throughput sequencing data. It has an automated installation script that sets up the code and third party tools used during analysis, and we’ve been working on improving the process to make getting started with […]','Motivation bcbio-nextgen is a community developed, best-practice pipeline for genomic data processing, performing automated variant calling and RNA-seq analyses from high throughput sequencing data. It has an automated installation script that sets up the code and third party tools used during analysis, and we’ve been working on improving the process to make getting started with […]','2014-03-06 00:00:00','2014-05-09 14:47:10.532000',0,'http://feedproxy.google.com/~r/bcbio/~3/e6pn3H4dKBk/');
INSERT INTO "planet_blogpost" VALUES(14,5,'http://bcbio.wordpress.com/?p=540','Updated comparison of variant detection methods: Ensemble, FreeBayes and minimal BAM preparation pipelines','Variant evaluation overview I previously discussed our approach for evaluating variant detection methods using a highly confident set of reference calls provided by NIST’s Genome in a Bottle consortium for the NA12878 human HapMap genome, In this post, I’ll update those conclusions based on recent improvements in GATK and FreeBayes. The comparisons use bcbio-nextgen, an […]','Variant evaluation overview I previously discussed our approach for evaluating variant detection methods using a highly confident set of reference calls provided by NIST’s Genome in a Bottle consortium for the NA12878 human HapMap genome, In this post, I’ll update those conclusions based on recent improvements in GATK and FreeBayes. The comparisons use bcbio-nextgen, an […]','2013-10-21 00:00:00','2014-05-09 14:47:10.540000',0,'http://feedproxy.google.com/~r/bcbio/~3/GItde4TSDGQ/');
INSERT INTO "planet_blogpost" VALUES(15,5,'http://bcbio.wordpress.com/?p=524','Summary from Bioinformatics Open Science Codefest 2013: Tools, infrastructure, standards and visualization','The 2013 Bioinformatics Open Source Conference (BOSC) starts tomorrow in Berlin, Germany. It’s a yearly conference devoted to community-based software development projects supporting biological research. Members of the Open Bioinformatics Foundation discuss implementations and approaches to better provide interoperable and reusable software, libraries and pipelines. For the past five years, a two day Codefest and […]','The 2013 Bioinformatics Open Source Conference (BOSC) starts tomorrow in Berlin, Germany. It’s a yearly conference devoted to community-based software development projects supporting biological research. Members of the Open Bioinformatics Foundation discuss implementations and approaches to better provide interoperable and reusable software, libraries and pipelines. For the past five years, a two day Codefest and […]','2013-07-18 00:00:00','2014-05-09 14:47:10.547000',0,'http://feedproxy.google.com/~r/bcbio/~3/9HApjLRjbrs/');
CREATE TABLE "planet_blog" ("feed" varchar(200) NOT NULL, "title" varchar(255) NOT NULL, "list_order" integer NOT NULL, "link" varchar(200) NOT NULL, "active" bool NOT NULL, "id" integer PRIMARY KEY, "desc" text NOT NULL);
INSERT INTO "planet_blog" VALUES('http://feeds.feedburner.com/MyWeblogOnBioinformaticsGenomeScienceNextGenerationSequencing','Kevin''s GATTACA World',0,'http://kevin-gattaca.blogspot.com/',1,1,'My Weblog on Bioinformatics, Genome Science and Next Generation Sequencing');
INSERT INTO "planet_blog" VALUES('http://feeds.feedburner.com/OmicsOmics','Omics! Omics!',0,'http://omicsomics.blogspot.com/',1,2,'A computational biologist''s personal views on new technologies &amp; publications on genomics &amp; proteomics and their impact on drug discovery');
INSERT INTO "planet_blog" VALUES('http://plindenbaum.blogspot.com/feeds/posts/default','YOKOFAKUN',0,'http://plindenbaum.blogspot.com/',1,3,'');
INSERT INTO "planet_blog" VALUES('http://nsaunders.wordpress.com/feed/','What You''re Doing Is Rather Desperate',0,'http://nsaunders.wordpress.com',1,4,'Notes from the life of a computational biologist');
INSERT INTO "planet_blog" VALUES('http://feeds2.feedburner.com/bcbio','Blue Collar Bioinformatics',0,'http://bcbio.wordpress.com',1,5,'');
CREATE TABLE "socialaccount_socialaccount" ("user_id" integer NOT NULL, "uid" varchar(255) NOT NULL, "last_login" datetime NOT NULL, "provider" varchar(30) NOT NULL, "extra_data" text NOT NULL, "id" integer PRIMARY KEY, "date_joined" datetime NOT NULL);
CREATE TABLE "socialaccount_socialapp_sites" ("id" integer NOT NULL PRIMARY KEY, "socialapp_id" integer NOT NULL, "site_id" integer NOT NULL);
INSERT INTO "socialaccount_socialapp_sites" VALUES(1,1,1);
INSERT INTO "socialaccount_socialapp_sites" VALUES(2,2,1);
INSERT INTO "socialaccount_socialapp_sites" VALUES(3,3,1);
CREATE TABLE "socialaccount_socialapp" ("key" varchar(100) NOT NULL, "secret" varchar(100) NOT NULL, "client_id" varchar(100) NOT NULL, "provider" varchar(30) NOT NULL, "id" integer PRIMARY KEY, "name" varchar(40) NOT NULL);
INSERT INTO "socialaccount_socialapp" VALUES('','','','persona',1,'persona');
INSERT INTO "socialaccount_socialapp" VALUES('','foobar','foobar','google',2,'google');
INSERT INTO "socialaccount_socialapp" VALUES('','foobar','foobar','github',3,'github');
CREATE TABLE "socialaccount_socialtoken" ("account_id" integer NOT NULL, "app_id" integer NOT NULL, "expires_at" datetime, "token" text NOT NULL, "id" integer PRIMARY KEY, "token_secret" text NOT NULL);
CREATE TABLE "djcelery_intervalschedule" ("id" integer NOT NULL PRIMARY KEY, "every" integer NOT NULL, "period" varchar(24) NOT NULL);
CREATE TABLE "djcelery_periodictasks" ("ident" smallint NOT NULL PRIMARY KEY, "last_update" datetime NOT NULL);
CREATE TABLE "djcelery_workerstate" ("id" integer NOT NULL PRIMARY KEY, "hostname" varchar(255) NOT NULL UNIQUE, "last_heartbeat" datetime NULL);
CREATE TABLE "djcelery_taskstate" ("id" integer NOT NULL PRIMARY KEY, "state" varchar(64) NOT NULL, "task_id" varchar(36) NOT NULL UNIQUE, "name" varchar(200) NULL, "tstamp" datetime NOT NULL, "args" text NULL, "kwargs" text NULL, "eta" datetime NULL, "expires" datetime NULL, "result" text NULL, "traceback" text NULL, "runtime" real NULL, "retries" integer NOT NULL, "worker_id" integer NULL, "hidden" bool NOT NULL);
CREATE TABLE "djcelery_periodictask" ("crontab_id" integer, "task" varchar(200) NOT NULL, "name" varchar(200) NOT NULL UNIQUE, "exchange" varchar(200), "args" text NOT NULL, "enabled" bool NOT NULL, "routing_key" varchar(200), "interval_id" integer, "last_run_at" datetime, "queue" varchar(200), "total_run_count" integer unsigned NOT NULL, "expires" datetime, "kwargs" text NOT NULL, "date_changed" datetime NOT NULL, "id" integer PRIMARY KEY, "description" text NOT NULL);
CREATE TABLE "celery_tasksetmeta" ("taskset_id" varchar(255) NOT NULL UNIQUE, "hidden" bool NOT NULL, "id" integer PRIMARY KEY, "date_done" datetime NOT NULL, "result" text NOT NULL);
CREATE TABLE "djcelery_crontabschedule" ("hour" varchar(64) NOT NULL, "day_of_month" varchar(64) NOT NULL, "day_of_week" varchar(64) NOT NULL, "month_of_year" varchar(64) NOT NULL, "id" integer PRIMARY KEY, "minute" varchar(64) NOT NULL);
CREATE TABLE "celery_taskmeta" ("status" varchar(50) NOT NULL, "task_id" varchar(255) NOT NULL UNIQUE, "date_done" datetime NOT NULL, "traceback" text, "meta" text NULL, "result" text, "hidden" bool NOT NULL, "id" integer PRIMARY KEY);
CREATE TABLE "djkombu_queue" ("id" integer NOT NULL PRIMARY KEY, "name" varchar(200) NOT NULL UNIQUE);
CREATE TABLE "djkombu_message" ("id" integer NOT NULL PRIMARY KEY, "visible" bool NOT NULL, "sent_at" datetime NULL, "payload" text NOT NULL, "queue_id" integer NOT NULL);
CREATE INDEX "auth_permission_37ef4eb4" ON "auth_permission" ("content_type_id");
CREATE INDEX "auth_group_permissions_5f412f9a" ON "auth_group_permissions" ("group_id");
CREATE INDEX "auth_group_permissions_83d7f98b" ON "auth_group_permissions" ("permission_id");
CREATE INDEX "messages_messagebody_e969df21" ON "messages_messagebody" ("author_id");
CREATE INDEX "messages_messagebody_cd3bbc30" ON "messages_messagebody" ("parent_msg_id");
CREATE INDEX "messages_message_6340c63c" ON "messages_message" ("user_id");
CREATE INDEX "messages_message_44a120f9" ON "messages_message" ("body_id");
CREATE INDEX "messages_message_403d8ff3" ON "messages_message" ("type");
CREATE INDEX "messages_message_bc4c5ddc" ON "messages_message" ("sent_at");
CREATE INDEX "django_admin_log_6340c63c" ON "django_admin_log" ("user_id");
CREATE INDEX "django_admin_log_37ef4eb4" ON "django_admin_log" ("content_type_id");
CREATE INDEX "django_flatpage_sites_872c4601" ON "django_flatpage_sites" ("flatpage_id");
CREATE INDEX "django_flatpage_sites_99732b5c" ON "django_flatpage_sites" ("site_id");
CREATE INDEX "django_flatpage_c379dc61" ON "django_flatpage" ("url");
CREATE INDEX "django_session_b7b81f0c" ON "django_session" ("expire_date");
CREATE INDEX "account_emailaddress_6340c63c" ON "account_emailaddress" ("user_id");
CREATE INDEX "account_emailconfirmation_a659cab3" ON "account_emailconfirmation" ("email_address_id");
CREATE INDEX "users_user_99732b5c" ON "users_user"("site_id");
CREATE UNIQUE INDEX "users_profile_tags_profile_id__tag_id" ON "users_profile_tags"("profile_id", "tag_id");
CREATE INDEX "users_profile_tags_5948a8a3" ON "users_profile_tags" ("profile_id");
CREATE INDEX "users_profile_tags_5659cca2" ON "users_profile_tags" ("tag_id");
CREATE UNIQUE INDEX "posts_post_tag_set_post_id__tag_id" ON "posts_post_tag_set"("post_id", "tag_id");
CREATE UNIQUE INDEX "posts_subscription_user_id__post_id" ON "posts_subscription"("user_id", "post_id");
CREATE INDEX "posts_tag_4da47e07" ON "posts_tag" ("name");
CREATE INDEX "posts_post_tag_set_87a49a9a" ON "posts_post_tag_set" ("post_id");
CREATE INDEX "posts_post_tag_set_5659cca2" ON "posts_post_tag_set" ("tag_id");
CREATE INDEX "posts_postview_87a49a9a" ON "posts_postview" ("post_id");
CREATE INDEX "posts_vote_e969df21" ON "posts_vote" ("author_id");
CREATE INDEX "posts_vote_87a49a9a" ON "posts_vote" ("post_id");
CREATE INDEX "posts_vote_403d8ff3" ON "posts_vote" ("type");
CREATE INDEX "posts_vote_eeede814" ON "posts_vote" ("date");
CREATE INDEX "posts_subscription_6340c63c" ON "posts_subscription" ("user_id");
CREATE INDEX "posts_subscription_87a49a9a" ON "posts_subscription" ("post_id");
CREATE INDEX "posts_subscription_403d8ff3" ON "posts_subscription" ("type");
CREATE INDEX "posts_subscription_eeede814" ON "posts_subscription" ("date");
CREATE INDEX "posts_post_e969df21" ON "posts_post"("author_id");
CREATE INDEX "posts_post_ce327ac6" ON "posts_post"("lastedit_user_id");
CREATE INDEX "posts_post_403d8ff3" ON "posts_post"("type");
CREATE INDEX "posts_post_b106351b" ON "posts_post"("vote_count");
CREATE INDEX "posts_post_8ffb2a07" ON "posts_post"("thread_score");
CREATE INDEX "posts_post_9f8d4624" ON "posts_post"("creation_date");
CREATE INDEX "posts_post_65fc64f1" ON "posts_post"("lastedit_date");
CREATE INDEX "posts_post_295a4710" ON "posts_post"("sticky");
CREATE INDEX "posts_post_82798ac5" ON "posts_post"("root_id");
CREATE INDEX "posts_post_410d0aac" ON "posts_post"("parent_id");
CREATE INDEX "posts_post_99732b5c" ON "posts_post"("site_id");
CREATE INDEX "posts_replytoken_6340c63c" ON "posts_replytoken" ("user_id");
CREATE INDEX "posts_replytoken_87a49a9a" ON "posts_replytoken" ("post_id");
CREATE INDEX "posts_emailentry_87a49a9a" ON "posts_emailentry" ("post_id");
CREATE INDEX "badges_award_6340c63c" ON "badges_award"("user_id");
CREATE INDEX "badges_award_540bf7e8" ON "badges_award"("badge_id");
CREATE INDEX "planet_blogpost_7894a4cc" ON "planet_blogpost" ("blog_id");
CREATE INDEX "planet_blogpost_9f8d4624" ON "planet_blogpost" ("creation_date");
CREATE INDEX "planet_blogpost_a889ce29" ON "planet_blogpost" ("insert_date");
CREATE INDEX "socialaccount_socialaccount_6340c63c" ON "socialaccount_socialaccount"("user_id");
CREATE UNIQUE INDEX "socialaccount_socialaccount_uid__provider" ON "socialaccount_socialaccount"("uid", "provider");
CREATE UNIQUE INDEX "socialaccount_socialapp_sites_socialapp_id__site_id" ON "socialaccount_socialapp_sites"("socialapp_id", "site_id");
CREATE INDEX "socialaccount_socialapp_sites_f2973cd1" ON "socialaccount_socialapp_sites" ("socialapp_id");
CREATE INDEX "socialaccount_socialapp_sites_99732b5c" ON "socialaccount_socialapp_sites" ("site_id");
CREATE UNIQUE INDEX "socialaccount_socialtoken_app_id__account_id" ON "socialaccount_socialtoken"("app_id", "account_id");
CREATE INDEX "djcelery_workerstate_11e400ef" ON "djcelery_workerstate" ("last_heartbeat");
CREATE INDEX "djcelery_taskstate_5654bf12" ON "djcelery_taskstate" ("state");
CREATE INDEX "djcelery_taskstate_4da47e07" ON "djcelery_taskstate" ("name");
CREATE INDEX "djcelery_taskstate_abaacd02" ON "djcelery_taskstate" ("tstamp");
CREATE INDEX "djcelery_taskstate_cac6a03d" ON "djcelery_taskstate" ("worker_id");
CREATE INDEX "djcelery_taskstate_2ff6b945" ON "djcelery_taskstate" ("hidden");
CREATE INDEX "djcelery_periodictask_7280124f" ON "djcelery_periodictask"("crontab_id");
CREATE INDEX "djcelery_periodictask_8905f60d" ON "djcelery_periodictask"("interval_id");
CREATE INDEX "djkombu_message_5907bb86" ON "djkombu_message" ("visible");
CREATE INDEX "djkombu_message_bc4c5ddc" ON "djkombu_message" ("sent_at");
CREATE INDEX "djkombu_message_c80a9385" ON "djkombu_message" ("queue_id");
COMMIT;
