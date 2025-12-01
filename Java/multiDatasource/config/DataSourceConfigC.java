
package com.example.multidatasource.config;

import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.boot.autoconfigure.orm.jpa.JpaProperties;
import org.springframework.boot.context.properties.ConfigurationProperties;
import org.springframework.boot.jdbc.DataSourceBuilder;
import org.springframework.boot.orm.jpa.EntityManagerFactoryBuilder;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.context.annotation.Primary;
import org.springframework.data.jpa.repository.config.EnableJpaRepositories;
import org.springframework.orm.jpa.JpaTransactionManager;
import org.springframework.orm.jpa.LocalContainerEntityManagerFactoryBean;
import org.springframework.transaction.PlatformTransactionManager;
import org.springframework.transaction.annotation.EnableTransactionManagement;

import javax.persistence.EntityManagerFactory;
import javax.sql.DataSource;
import java.util.HashMap;
import java.util.Map;

@Configuration
@EnableJpaRepositories(
    basePackages = "com.example.multidatasource.repository.c",
    entityManagerFactoryRef = "entityManagerFactoryC",
    transactionManagerRef = "transactionManagerC"
)
public class DataSourceConfigC {
    @Bean(name = "dataSourceC")
    @ConfigurationProperties(prefix = "spring.datasource.c")
    public DataSource dataSourceC() {
        return DataSourceBuilder.create().build();
    }

    @Bean(name = "jpaPropertiesC")
    @ConfigurationProperties(prefix = "spring.jpa.c")
    public JpaProperties jpaPropertiesC() {
        return new JpaProperties();
    }

    @Bean(name = "entityManagerFactoryC")
    public LocalContainerEntityManagerFactoryBean entityManagerFactoryC(
            EntityManagerFactoryBuilder builder,
            @Qualifier("dataSourceC") DataSource dataSource) {
        Map<String, String> props = jpaPropertiesB().getProperties();

        return builder
                .dataSource(dataSource)
                .packages("com.example.multidatasource.entity.c")
                .persistenceUnit("c")
                .properties(props())
                .build();
    }

    @Bean(name = "transactionManagerC")
    public PlatformTransactionManager transactionManagerC(
            @Qualifier("entityManagerFactoryC") EntityManagerFactory entityManagerFactory) {
        return new JpaTransactionManager(entityManagerFactory);
    }
}
